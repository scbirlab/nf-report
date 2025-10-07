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


process fetch_fasta_from_uniprot_id {

   tag "${uniprot_id}"

   publishDir( 
      "${params.outputs}/sequences", 
      mode: 'copy',
      saveAs: { "${id}-${uniprot_id}.fasta" },
   )

   input:
   tuple val( id ), val( uniprot_id )

   output:
   tuple val( id ), path( "protein.fasta" )

   script:
   """
   wget "https://rest.uniprot.org/uniprotkb/stream?query=(accession:${uniprot_id})&format=fasta&download=true&compressed=false" \
      -O protein.fasta \
   || (
      echo "Failed to download Uniprot ID ${uniprot_id} from UniProt"
      exit 1   
   )
   """

}


process fetch_fastas_from_uniprot_ids {

   tag "${id[0]}...${id[-1]}"

   errorStrategy 'retry'
   maxRetries 2

   publishDir( 
      "${params.outputs}/sequences", 
      mode: 'copy',
      saveAs: { "${id[0]}-${id[-1]}.fasta" },
   )

   input:
   val id

   output:
   path "proteins.fasta"

   script:
   """
   set -x
   curl -X GET --header 'Accept:text/x-fasta' \
      'https://www.ebi.ac.uk/proteins/api/proteins?offset=0&size=-1&accession=${id.join(',')}' \
   > proteins.fasta

   """

}

process fetch_species_gene_names {

   tag "${id}:${column}"

   publishDir( 
      "${params.outputs}/targets", 
      mode: 'copy',
      saveAs: { "${id}-${it}" },
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

   printf '${column}\\tspecies_target_name\\tspecies_target_locus\\n' \
   > targets0.tsv
   
   if ls ids_* 1> /dev/null 2>&1
   then
      url='https://www.ebi.ac.uk/proteins/api/proteins'
      base_query='offset=0&size=-1'
      header='Accept:application/json'

      
      for f in ids_*
      do
         these_ids=\$(tr \$'\\n' , < "\$f")
         curl -s -X GET --header \$header \
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


process map_uniprot_ids_from_file {

   tag "${id}-${column}:${from_type}"

   publishDir( 
      "${params.outputs}/uniprot_map", 
      mode: 'copy',
      saveAs: { "${id}-${column}.txt" },
   )

   input:
   tuple val( id ), val( organism_id ), path( filename ), val( from_type ), val( column )

   output:
   tuple val( id ), path( "uniprot-ids.txt" )

   script:
   """
   bash ${projectDir}/bin/uniprot/uniprot-ids.sh \
      "${filename}" "${column}" \
      "uniprot-ids0.txt" \
      UniProtKB ${from_type} \
      "${organism_id}"
   cut -f2 -d, < "uniprot-ids0.txt" \
   | tail -n+2 \
   | grep -v '^\$' \
   > "uniprot-ids.txt"
   """

}


process map_gene_names_from_file {

   tag "${id}-${column}"

   publishDir( 
      "${params.outputs}/ppi", 
      mode: 'copy',
      saveAs: { "${id}-${table.getSimpleName()}-${column}.txt" },
   )

   input:
   tuple val( id ), path( table )
   val column

   output:
   tuple val( id ), path( "named-ids.tsv" )

   script:
   """
   set -x

   get_column_number () (
      local col_name="\$1"
      head -n1 | tr \$'\t' \$'\n' | grep -n "\$col_name" | cut -d: -f1
   )

   sort_table () (
      local col_name="\$1"
      local tempfile="aaaaaa"
      cat > \$tempfile

      col_number=\$(get_column_number "\$col_name" < \$tempfile)
      head -n1 \$tempfile \
      | cat - <(tail -n+2 \$tempfile | sort -k"\$col_number") \
      && rm \$tempfile
   )

   bash ${projectDir}/bin/uniprot/uniprot-ids.sh \
      "${table}" "${column}" \
      "uniprot-ids0.txt" \
      Gene_Name UniProtKB
   sort_table "${column}" < "uniprot-ids0.txt" \
   > "uniprot-ids.txt"
   join_col=\$(get_column_number "${column}" < "${table}")
   join --header -t \$'\t' \
      -1 1 -2 "\$join_col"  \
      <(tr , \$'\t' < "uniprot-ids.csv") \
      <(sort_table "${column}" < "${table}") \
   > "named-ids.tsv"
   """

   stub:
   """

   set -x

   get_column_number () (
      local col_name="\$1"
      head -n1 | tr \$'\t' \$'\n' | grep -n "\$col_name" | cut -d: -f1
   )

   sort_table () (
      local col_name="\$1"
      local tempfile="aaaaaa"
      cat > \$tempfile

      col_number=\$(get_column_number "\$col_name" < \$tempfile)
      head -n1 \$tempfile \
      | cat - <(tail -n+2 \$tempfile | sort -k"\$col_number") \
      && rm \$tempfile
   )

   cp ${projectDir}/data/559292-dca-stub.tsv input.tsv

   bash ${projectDir}/bin/uniprot/uniprot-ids.sh \
      input.tsv "${column}" \
      "uniprot-ids0.txt" \
      Gene_Name UniProtKB
   sort_table "${column}" < "uniprot-ids0.txt" \
   | sed 's/Gene_Name/${column}_gene_name/g' \
   > "uniprot-ids.csv"
   join_col=\$(get_column_number "${column}" < "input.tsv")
   join --header -t \$'\t' \
      -1 1 -2 "\$join_col"  \
      <(tr , \$'\t' < "uniprot-ids.csv") \
      <(sort_table "${column}" < "input.tsv") \
   > "named-ids.tsv"

   """

}