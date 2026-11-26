process fetch_fastas_from_organism_id {

   errorStrategy { task.exitStatus == 35 ? 'retry' : ( task.exitStatus in [45, 46] ? 'ignore' : 'terminate') }  // sometimes UniProt fails to respond
   maxRetries 1
   stageInMode 'link'
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

   function get_proteome_id() {
      curl "https://rest.uniprot.org/proteomes/search?query=(taxonomy_id:${organism_id})&format=json" \
      | jq -r '
         .results[] 
         | select(.proteomeType == "'"\$1"' proteome") 
         | .id
      ' | head -n1
   }

   function has_fasta() {
      [ -s "\$1" ] && gunzip -t "\$1" && zcat "\$1" \
      | grep '^>' > /dev/null
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

   # Returns 0 = usable FASTA written, 1 = reached the server but nothing usable
   # (e.g. purged/404), 2 = connection or server-side trouble, worth a retry.
   function try_fetch() {
      local status
      status=\$(curl -s -o "\$2" -w '%{http_code}' "\$1" || echo "000")
      if has_fasta "\$2"
      then
         return 0
      fi
      case "\$status" in
         000|5??) return 2 ;;
         *)       return 1 ;;
      esac
   }

   QUERIES=("Reference and representative" "Reference" "Representative" "Non Reference" "Other")
   PROTEOME_ID=
   for q in "\${QUERIES[@]}"
   do
      PROTEOME_ID=\$(get_proteome_id "\$q")
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

   TROUBLE=0

   if try_fetch "https://rest.uniprot.org/uniprotkb/stream?query=(proteome:\$PROTEOME_ID)&format=fasta&download=true&compressed=true" uniprotkb.fasta.gz
   then
      mv uniprotkb.fasta.gz proteome.fasta.gz
      exit 0
   elif [ \$? -eq 2 ]
   then
      TROUBLE=1
   fi

   echo "Proteome \$PROTEOME_ID not recoverable from UniProtKB (taxonomy ${organism_id}); trying UniParc"

   if try_fetch "https://rest.uniprot.org/uniparc/stream?query=(upid:\$PROTEOME_ID)&format=fasta&compressed=true" uniparc.fasta.gz
   then
      normalize_headers uniparc.fasta.gz proteome.fasta.gz
      exit 0
   elif [ \$? -eq 2 ]
   then
      TROUBLE=1
   fi

   if [ "\$TROUBLE" -eq 1 ]
   then
      echo "UniProt/UniParc unreachable for proteome \$PROTEOME_ID (taxonomy ${organism_id})"
      exit 35
   else
      echo "Proteome \$PROTEOME_ID has no recoverable sequences in UniProtKB or UniParc (taxonomy ${organism_id})"
      exit 46
   fi

   """

}


process fetch_fastas_from_uniprot_ids {

   tag "${id[0]}...${id[-1]}"

   errorStrategy 'retry'
   maxRetries 2

   publishDir( 
      "${params.outputs}/sequences", 
      mode: 'copy',
      saveAs: { "${id[0]}.${id[-1]}.fasta.gz" },
   )

   input:
   val id

   output:
   path "proteins.fasta.gz"

   script:
   """
   set -x
   curl -X GET --header 'Accept:text/x-fasta'  -A 'scbirlab-nf-report/0.4 (+https://scbirlab.org; contact: eachan.johnson@crick.ac.uk)' \
      'https://www.ebi.ac.uk/proteins/api/proteins?offset=0&size=-1&accession=${id.join(',')}' \
   > proteins.fasta
   gzip --best proteins.fasta

   """

}

process fetch_species_gene_names {

   tag "${id}:${column}"
   stageInMode 'link'
   errorStrategy { if(task.attempt > 2) { return 'retry' } else { return 'ignore' } }  // retry if network issue, ignore if input issue
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
   set -x

   col=\$(zcat "${table}" | head -n1  | tr \$'\\t' \$'\\n' | grep -nFw "${column}" | cut -d: -f1)
   zcat "${table}" | tail -n+2 | cut -f"\$col" | sort -u | split -l50 - 'ids_'

   printf '${column}\\tortholog_target_name\\tortholog_target_locus\\n' \
   > targets0.tsv
   
   if ls ids_* 1> /dev/null 2>&1
   then
      url='https://www.ebi.ac.uk/proteins/api/proteins'
      base_query='offset=0&size=-1'
      header='Accept:application/json'

      
      for f in ids_*
      do
         these_ids=\$(tr \$'\\n' , < "\$f")
         curl -s -X GET --header \$header  -A 'scbirlab-nf-report/0.4 (+https://scbirlab.org; contact: eachan.johnson@crick.ac.uk)' \
            "\${url}?\${base_query}&accession=\${these_ids}" \
         | jq -r '
            .[] | [
               .accession, 
               (.gene[0].name.value // "NA"), 
               ((.gene[0].olnNames // (.gene[0].orfNames // []))[0].value // "NA")
            ] | @tsv
         ' \
         >> targets0.tsv
      done
   else
      echo "" > targets0.tsv
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


process fetch_species_gene_names2 {

   tag "${id}:${column}"
   stageInMode 'link'
   errorStrategy { if(task.attempt > 2) { return 'retry' } else { return 'ignore' } }  // retry if network issue, ignore if input issue
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
   set -x

   col=\$(zcat "${table}" | head -n1 | tr \$'\\t' \$'\\n' | grep -nFw "${column}" | cut -d: -f1)
   zcat "${table}" | tail -n+2 | cut -f"\$col" | sort -u > all_ids.txt

   # UniParc accessions (proteome archived off UniProtKB) need a different endpoint
   # and a taxon-specific filter; everything else is a UniProtKB accession.
   grep -E '^UPI[0-9A-Fa-f]{10}\$' all_ids.txt > uniparc_ids.txt || true
   grep -vE '^UPI[0-9A-Fa-f]{10}\$' all_ids.txt > uniprotkb_ids.txt || true
   split -l50 uniprotkb_ids.txt 'ids_'

   printf '${column}\\tortholog_target_name\\tortholog_target_locus\\tortholog_name_source_db\\tortholog_name_source_id\\tgene_name_source\\n' \
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
         curl -s -X GET --header "\$header" -A "\$ua" \
            "\${url}?\${base_query}&accession=\${these_ids}" \
         | jq -r '
            .[] | [
            .accession, 
            (.gene[0].name.value // "NA"), 
            ((.gene[0].olnNames // (.gene[0].orfNames // []))[0].value // "NA"),
            "NA",
            .accession,
            "uniprotkb"
         ] | @tsv
         ' \
         >> targets0.tsv
      done
   fi

   # --- UniParc accessions ---
   if [ -s uniparc_ids.txt ]
   then
      while read -r upi
      do
         resp=\$(curl -s -X GET --header \$header -A "\$ua" \
            "\${uniparc_url}/\${upi}?rfTaxId=${id}")
         if [ -z "\$resp" ]
         then
            printf '%s\\tNA\\tNA\\tNA\\tNA\\tuniparc_lookup_failed\\n' "\$upi" \
            >> targets0.tsv
            continue
         fi
         echo "\$resp" \
         | jq -r \
            --arg acc "\$upi" '
            (.dbReference // []) as \$refs
            | ( \$refs | map(select(.type == "UniProtKB/TrEMBL")) | first ) as \$trembl
            | ( \$trembl // \$refs[0] // {} ) as \$entry
            | ( \$entry.property // [] ) as \$props
            | ( [\$props[] | select(.type=="gene_name") | .value] | first // "NA" ) as \$gene
            | ( [\$props[] | select(.type=="protein_name") | .value] | first // "NA" ) as \$protein
            | [
                \$acc, 
                \$gene, 
                \$protein,
                (\$entry.type // "NA"), 
                (\$entry.id // "NA"),
                (if (\$refs | length) == 0 then "uniparc_no_taxon_match" else "uniparc" end)
              ] | @tsv
         ' >> targets0.tsv
      done < uniparc_ids.txt
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
