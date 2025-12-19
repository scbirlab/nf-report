process fetch_fastas_from_organism_id {

   errorStrategy { task.exitStatus == 35 ? 'retry' : ( task.exitStatus == 45 ? 'ignore' : 'terminate') }  // sometimes UniProt fails to respond
   maxRetries 1
   stageInMode 'link'

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
   QUERIES=("Reference and representative" "Reference" "Representative" "Other")
   PROTEOME_ID=
   for q in "\${QUERIES[@]}"
   do
      PROTEOME_ID=\$(get_proteome_id "\$q")
      if [ -n "\$PROTEOME_ID" ]
      then 
         break
      fi
   done

   if [ -n "\$PROTEOME_ID" ]
   then
      curl -v "https://rest.uniprot.org/uniprotkb/stream?query=(proteome:\$PROTEOME_ID)&format=fasta&download=true&compressed=true" \
      > proteome.fasta.gz
   else
      echo "Failed to download taxonomy ID ${organism_id} with proteome ID \$PROTEOME_ID from UniProt"
      exit 45
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
