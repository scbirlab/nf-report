process fetch_fastas_from_organism_id {

   errorStrategy 'retry'  // sometimes UniProt fails to respond
   maxRetries 2

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
   function get_proteome_id() {
      curl -s "https://rest.uniprot.org/proteomes/search?query=(taxonomy_id:${organism_id})&format=json" \
      | jq -r '.results[] | select(.proteomeType == "'"\$1"' proteome").id'
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
      wget "https://rest.uniprot.org/uniprotkb/stream?query=(proteome:\$PROTEOME_ID)&format=fasta&download=true&compressed=true" \
      -O proteome.fasta.gz
   else
      echo "Failed to download taxonomy ID ${organism_id} with proteome ID \$PROTEOME_ID from UniProt"
      exit 1
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
      saveAs: { "${id[0]}.${id[-1]}.fasta" },
   )

   input:
   val id

   output:
   path "proteins.fasta"

   script:
   """
   set -x
   curl -X GET --header 'Accept:text/x-fasta'  -A 'scbirlab-nf-report/0.4 (+https://scbirlab.org; contact: eachan.johnson@crick.ac.uk)' \
      'https://www.ebi.ac.uk/proteins/api/proteins?offset=0&size=-1&accession=${id.join(',')}' \
   > proteins.fasta

   """

}

process fetch_species_gene_names {

   tag "${id}:${column}"

   publishDir( 
      "${params.outputs}/targets", 
      mode: 'copy',
      saveAs: { "${id}.${it}" },
   )

   input:
   tuple val( id ), path( table ), path( taxon_table )
   val column

   output:
   tuple val( id ), path( 'targets.tsv' )

   script:
   """
   set -x

   col=\$(head -n1 "${table}" | tr \$'\\t' \$'\\n' | grep -nFw "${column}" | cut -d: -f1)
   tail -n+2 "${table}" | cut -f"\$col" | sort -u | split -l50 - 'ids_'

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
   
   (
      pd.read_csv("${taxon_table}", sep="\\t")
      .merge(
         pd.read_csv("${table}", sep="\\t"),
      )
      .merge(
         pd.read_csv("targets0.tsv", sep="\\t"),
      )
      .drop_duplicates()
      .to_csv("targets.tsv", sep="\\t", index=False)
   )
   
   '

   """

}
