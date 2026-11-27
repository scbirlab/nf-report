process fetch_fastas_from_organism_id {

   errorStrategy { task.exitStatus == 35 ? 'retry' : ( task.exitStatus in [45, 46] ? 'ignore' : 'terminate') }  // sometimes UniProt fails to respond or rate-limits
   maxRetries 3
   time '1h'

   tag "${id}"

   publishDir( 
      "${params.outputs}/proteome-sequences", 
      mode: 'copy',
      saveAs: { "${organism_id}.fasta.gz" }
   )

   input:
   tuple val( id ), val( organism_id )

   output:
   tuple val( id ), path( "proteome.fasta.gz" )

   script:
   """
   set -eox pipefail

   UA='scbirlab-nf-report/0.4 (+https://scbirlab.org)'

   # Back off on retries (UniProt rate-limits with HTTP 429)
   sleep \$(( (${task.attempt} - 1) * 30 ))

   # Prints the first proteome ID of the given type, or nothing if there is none.
   # Returns 2 if UniProt could not be reached or returned an error.
   function get_proteome_id() {
      # -o (not stdout) so that a retry rewrites the file instead of appending to a partial body
      curl -sS --fail --retry 3 --retry-delay 5 -A "\$UA" -o proteome-search.json \
         "https://rest.uniprot.org/proteomes/search?query=(taxonomy_id:${organism_id})&format=json" \
      || return 2
      jq -r --arg type "\$1 proteome" '
         first(.results[] | select(.proteomeType == \$type) | .id) // empty
      ' proteome-search.json
   }

   function normalize_headers() {
      zcat "\$1" \
      | awk '
         /^>/ {
            header = substr(\$0, 2)
            n = split(header, parts, "|")
            acc = (n >= 2) ? parts[2] : header
            sub(/ .*/, "", acc)
            print ">ref|" acc "|" acc
            next
         }
         { print \$0 }
      ' \
      | gzip -c > "\$2"
   }

   # Returns
   #   0 = usable FASTA written
   #   1 = confirmed absence (HTTP 200 with nothing in it, or 404/410)
   #   2 = transient: connection failure, rate-limit, server error, or a 200 whose
   #       body is not a complete gzip (truncated/corrupt download)
   #   3 = unexpected permanent error (any other HTTP status)
   # Status is checked first, so a bad response can never be mistaken for absence.
   function try_fetch() {
      local status
      # on connection failure curl exits non-zero but -w has already printed 000
      status=\$(curl -s -A "\$UA" -o "\$2" -w '%{http_code}' "\$1") || status=000
      case "\$status" in
         200)
            [ -s "\$2" ] || return 1
            gunzip -t "\$2" 2> /dev/null || return 2
            zcat "\$2" | grep '^>' > /dev/null && return 0
            return 1
            ;;
         404|410)     return 1 ;;
         000|429|5??) return 2 ;;
         *)           echo "Unexpected HTTP \$status from \$1" >&2; return 3 ;;
      esac
   }

   QUERIES=("Reference and representative" "Reference" "Representative" "Non Reference" "Other")
   PROTEOME_ID=
   for q in "\${QUERIES[@]}"
   do
      rc=0
      PROTEOME_ID=\$(get_proteome_id "\$q") || rc=\$?
      if [ "\$rc" -ne 0 ]
      then
         echo "UniProt proteome search failed for taxonomy ID ${organism_id}"
         exit 35
      fi
      if [ -n "\$PROTEOME_ID" ]
      then 
         echo "Found \$PROTEOME_ID"
         break
      fi
   done

   if [ -z "\$PROTEOME_ID" ]
   then
      echo "No proteome catalogued at all for taxonomy ID ${organism_id}"
      exit 45
   fi

   # A transient failure (rc 2) must retry the whole task, not fall through to
   # UniParc: that would silently swap UniProtKB accessions for UPI accessions.
   rc=0
   try_fetch "https://rest.uniprot.org/uniprotkb/stream?query=(proteome:\$PROTEOME_ID)&format=fasta&download=true&compressed=true" uniprotkb.fasta.gz || rc=\$?
   if [ "\$rc" -eq 0 ]
   then
      mv uniprotkb.fasta.gz proteome.fasta.gz
      exit 0
   elif [ "\$rc" -eq 2 ]
   then
      echo "UniProtKB unreachable, rate-limiting or returned a corrupt download for proteome \$PROTEOME_ID (taxonomy ${organism_id})"
      exit 35
   elif [ "\$rc" -eq 3 ]
   then
      echo "Unexpected UniProtKB response for proteome \$PROTEOME_ID (taxonomy ${organism_id})"
      exit 1
   fi

   echo "Proteome \$PROTEOME_ID not recoverable from UniProtKB (taxonomy ${organism_id}); trying UniParc"

   rc=0
   try_fetch "https://rest.uniprot.org/uniparc/stream?query=(upid:\$PROTEOME_ID)&format=fasta&compressed=true" uniparc.fasta.gz || rc=\$?
   if [ "\$rc" -eq 0 ]
   then
      normalize_headers uniparc.fasta.gz proteome.fasta.gz
      exit 0
   elif [ "\$rc" -eq 2 ]
   then
      echo "UniParc unreachable, rate-limiting or returned a corrupt download for proteome \$PROTEOME_ID (taxonomy ${organism_id})"
      exit 35
   elif [ "\$rc" -eq 3 ]
   then
      echo "Unexpected UniParc response for proteome \$PROTEOME_ID (taxonomy ${organism_id})"
      exit 1
   fi

   echo "Proteome \$PROTEOME_ID has no recoverable sequences in UniProtKB or UniParc (taxonomy ${organism_id})"
   exit 46

   """

}


process fetch_fastas_from_uniprot_ids {

   tag "${id[0]}...${id[-1]}"

   errorStrategy 'retry'
   maxRetries 2

   publishDir( 
      "${params.outputs}/sequences", 
      mode: 'copy',
      saveAs: { "${id[0]}.${id[-1]}.${it}" },
   )

   input:
   val id

   output:
   path "fasta.gz", emit: fasta
   path "missing-accessions.txt", emit: missing
   path "unexpected-accessions.txt", emit: unexpected

   script:
   """
   set -euxo pipefail

   cat > requested.txt <<'EOF_IDS'
${id.join('\n')}
EOF_IDS
   sort -u requested.txt > requested.sorted.txt

   # --fail makes HTTP errors a non-zero exit, so errorStrategy 'retry' actually applies
   curl -sS --fail --retry 3 --retry-delay 5 --header 'Accept:text/x-fasta' -A 'scbirlab-nf-report/0.4 (+https://scbirlab.org; contact: eachan.johnson@crick.ac.uk)' \
      -o proteins.fasta \
      'https://www.ebi.ac.uk/proteins/api/proteins?offset=0&size=-1&accession=${id.join(',')}'

   # The API silently omits accessions it does not know (obsolete, demerged, isoforms),
   # so compare what came back against what was asked for.
   { grep '^>' proteins.fasta || true; } | cut -d'|' -f2 | sort -u > returned.txt
   comm -23 requested.sorted.txt returned.txt > missing-accessions.txt
   comm -13 requested.sorted.txt returned.txt > unexpected-accessions.txt

   # Structural checks: duplicated records (e.g. a retry appended to a partial body) or
   # headers with no sequence mean the download is corrupt, not merely incomplete.
   { grep '^>' proteins.fasta || true; } | cut -d'|' -f2 | sort | uniq -d > duplicated.txt
   if [ -s duplicated.txt ]
   then
      echo "ERROR: duplicated FASTA records for \$(wc -l < duplicated.txt) accessions" >&2
      exit 1
   fi
   empty_records=\$(awk '/^>/ { if (hdr && !seq) n++; hdr=1; seq=0; next } NF { seq=1 } END { if (hdr && !seq) n++; print n+0 }' proteins.fasta)
   if [ "\$empty_records" -gt 0 ]
   then
      echo "ERROR: \$empty_records FASTA headers have no sequence" >&2
      exit 1
   fi

   if [ ! -s returned.txt ]
   then
      echo "ERROR: no sequences returned for a batch of \$(wc -l < requested.sorted.txt) accessions" >&2
      exit 1
   fi
   if [ -s unexpected-accessions.txt ]
   then
      # e.g. secondary accessions resolved to a primary one; informative unless strict
      echo "WARNING: \$(wc -l < unexpected-accessions.txt) returned accessions were not requested; see unexpected-accessions.txt" >&2
      if [ "${params.fail_on_missing_sequences}" == "true" ]
      then
         exit 1
      fi
   fi
   if [ -s missing-accessions.txt ]
   then
      echo "WARNING: \$(wc -l < missing-accessions.txt) of \$(wc -l < requested.sorted.txt) accessions not returned; see missing-accessions.txt" >&2
      if [ "${params.fail_on_missing_sequences}" == "true" ]
      then
         exit 1
      fi
   fi
   gzip --best -c proteins.fasta > fasta.gz

   """

}

process fetch_species_gene_names {

   tag "${id}:${column}"
   stageInMode 'link'
   errorStrategy { task.exitStatus == 35 && task.attempt <= task.maxRetries ? 'retry' : 'terminate' }  // 35 = transient network failure (retry); anything else, or retries exhausted, stops the run rather than silently dropping an organism
   maxRetries 2

   publishDir( 
      "${params.outputs}/targets/tables", 
      mode: 'copy',
      saveAs: { "${id}.${it}" },
   )

   input:
   tuple val( id ), path( table ), path( taxon_table )
   val column

   output:
   tuple val( id ), path( 'targets.tsv.gz' )

   script:
   """
   set -euxo pipefail

   # Exact (whole-header) match for the accession column
   IFS= read -r table_header < <(zcat "${table}")
   col=\$(
      awk -F'\\t' -v c="${column}" '{ for ( i=1; i<=NF; i++ ) if ( \$i == c ) { print i; exit } }' \
      <<< "\$table_header"
   )
   if [ -z "\$col" ]
   then
      echo "Column ${column} not found in ${table}" >&2
      exit 1
   fi
   zcat "${table}" | tail -n+2 | cut -f"\$col" | sort -u | awk 'NF' > all_ids.txt

   # UniParc accessions (proteome archived off UniProtKB) need a different endpoint
   # and a taxon-specific filter; everything else is a UniProtKB accession.
   grep -E '^UPI[0-9A-Fa-f]{10}\$' all_ids.txt > uniparc_ids.txt || true
   grep -vE '^UPI[0-9A-Fa-f]{10}\$' all_ids.txt > uniprotkb_ids.txt || true
   split -l50 uniprotkb_ids.txt 'ids_'

   printf '${column}\\tortholog_target_name\\tortholog_target_locus\\tortholog_target_protein_name\\tortholog_name_source_db\\tortholog_name_source_id\\tgene_name_source\\n' \
   > targets0.tsv

   url='https://www.ebi.ac.uk/proteins/api/proteins'
   uniparc_url='https://www.ebi.ac.uk/proteins/api/uniparc/upi'
   base_query='offset=0&size=-1'
   header='Accept:application/json'
   ua='scbirlab-nf-report/0.4 (+https://scbirlab.org; contact: eachan.johnson@crick.ac.uk)'

   # --- UniProtKB accessions ---
   if ls ids_* 1> /dev/null 2>&1
   then
      for f in ids_*
      do
         these_ids=\$(tr \$'\\n' , < "\$f")
         status=\$(curl -sS --retry 3 --retry-delay 5 --header "\$header" -A "\$ua" \
            -o response.json -w '%{http_code}' \
            "\${url}?\${base_query}&accession=\${these_ids}") || status=000
         case "\$status" in
            200) ;;
            000|429|5??) echo "Transient HTTP \$status from UniProt Proteins API" >&2; exit 35 ;;
            *) echo "Permanent HTTP \$status from UniProt Proteins API for batch \$f" >&2; exit 1 ;;
         esac
         jq -e 'type == "array"' response.json > /dev/null || exit 35
         jq -r '
            ( if type == "array" then .[] else empty end )
            | [
               .accession, 
               (.gene[0].name.value // "NA"), 
               ((.gene[0].olnNames // (.gene[0].orfNames // []))[0].value // "NA"),
               (.protein.recommendedName.fullName.value // .protein.submittedName[0].fullName.value // "NA"),
               "NA",
               .accession,
               "uniprotkb"
            ] | @tsv
         ' response.json \
         >> targets0.tsv
      done
   fi

   # --- UniParc accessions (one request per UPI, a few at a time) ---
   if [ -s uniparc_ids.txt ]
   then
      mkdir uniparc_rows

      lookup_upi() {
         local upi="\$1"
         local out="uniparc_rows/\${upi}.tsv"
         local resp="uniparc_rows/\${upi}.json"
         local status
         status=\$(curl -sS --retry 3 --retry-delay 5 --header "\$header" -A "\$ua" \
            -o "\$resp" -w '%{http_code}' "\${uniparc_url}/\${upi}?rfTaxId=${id}") || status=000
         case "\$status" in
            200) ;;
            000|429|5??)
               # exhausted transient failure: same policy as UniProtKB, stop rather than drop annotations
               echo "UniParc lookup of \$upi failed transiently (HTTP \$status)" | tee "uniparc_rows/\${upi}.transient" >&2
               return 0
               ;;
            404|410)
               # the server answered but has nothing for this UPI/taxon
               printf '%s\\tNA\\tNA\\tNA\\tNA\\tNA\\tuniparc_lookup_failed\\n' "\$upi" > "\$out"
               return 0
               ;;
            *)
               # 400/401/403/...: a malformed or refused query would fail for every UPI, so stop
               # (xargs then exits 123 and set -e ends the task)
               echo "Unexpected HTTP \$status from UniParc for \$upi" >&2
               return 1
               ;;
         esac
         if ! jq -e 'type == "object"' "\$resp" > /dev/null 2>&1
         then
            # a 200 that is not the expected JSON object is a bad download: treat as transient
            echo "UniParc lookup of \$upi returned an unexpected body" | tee "uniparc_rows/\${upi}.transient" >&2
            return 0
         fi
         jq -r --arg acc "\$upi" '
            (.dbReference // []) as \$refs
            | ( \$refs | map(select(.type == "UniProtKB/TrEMBL")) | first ) as \$trembl
            | ( \$trembl // \$refs[0] // {} ) as \$entry
            | ( \$entry.property // [] ) as \$props
            | ( [\$props[] | select(.type=="gene_name") | .value] | first // "NA" ) as \$gene
            | ( [\$props[] | select(.type=="protein_name") | .value] | first // "NA" ) as \$protein
            | [
                \$acc, 
                \$gene, 
                "NA",
                \$protein,
                (\$entry.type // "NA"), 
                (\$entry.id // "NA"),
                (if (\$refs | length) == 0 then "uniparc_no_taxon_match" else "uniparc" end)
              ] | @tsv
         ' "\$resp" > "\$out"
      }
      export -f lookup_upi
      export header ua uniparc_url

      # xargs -P 0 means unlimited concurrency, so insist on a positive integer
      if ! [[ "${params.uniparc_workers}" =~ ^[1-9][0-9]*\$ ]]
      then
         echo "params.uniparc_workers must be a positive integer, got '${params.uniparc_workers}'" >&2
         exit 1
      fi
      # strict mode in the worker shells too; a worker failure makes xargs exit 123 and ends the task
      xargs -a uniparc_ids.txt -P ${params.uniparc_workers} -I{} bash -euo pipefail -c 'lookup_upi "\$1"' _ {}
      if ls uniparc_rows/*.transient 1> /dev/null 2>&1
      then
         echo "UniParc lookups failed transiently for \$(ls uniparc_rows/*.transient | wc -l) accessions:" >&2
         cat uniparc_rows/*.transient >&2
         exit 35
      fi
      cat uniparc_rows/*.tsv >> targets0.tsv
   fi

   python -c '
   import pandas as pd
   import numpy as np
   (
      pd.read_csv("${taxon_table}", sep="\\t")
      .merge(
         pd.read_csv("${table}", sep="\\t"),
      )
      .merge(
         pd.read_csv("targets0.tsv", sep="\\t"),
         how="left",
      )
      .drop_duplicates()
      .assign(
         ortholog_taxon_id="${id}",
         target_is_human=lambda x: x["target_taxon_id"] == 9606, 
         target_is_bacteria=lambda x: x["target_taxon_l1"] == "Bacteria",
         ortholog_target_name=lambda x: np.where(
               x["ortholog_target_name"].isna(), 
               x["ortholog_target_locus"], 
               x["ortholog_target_name"],
         ),
      )
      .to_csv("targets.tsv.gz", sep="\\t", index=False)
   )
   '

   """

}
