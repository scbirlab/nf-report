process make_diamond_db {

    errorStrategy 'ignore'

    tag "${id}"

    input:
    tuple val( id ), path( fasta )

    output:
    tuple val( id ), path( "db.dmnd" )

    script:
    """
    if [ "${fasta}" == "*.gz" ]
    then
        FILE_SIZE=\$(gzip -l "${fasta}" | awk 'NR==2 {print \$2}')
    elif [ -s "${fasta}" ]
    then
        FILE_SIZE=1
    else
        FILE_SIZE=0
    fi

    if [ "\$FILE_SIZE" -gt 0 ]
    then
        diamond makedb --in "${fasta}" -d db.dmnd
    else     
        echo "${fasta} is empty!"
        exit 1
    fi
    """

}


process diamond_blastp {

    tag "${id}"

    publishDir( 
        "${params.outputs}/blast-results", 
        mode: 'copy',
        saveAs: { "${id}-${it}" }
    )

    input:
    tuple val( id ), path( db ), path( queries )

    output:
    tuple val( id ), path( "hits.tsv" ), emit: data
    tuple val( id ), path( "hit_count.tsv" ), emit: stats

    script:
    """
    header="target_uniprot_id,ortholog_uniprot_id,target_accession,ortholog_accession,target_length,ortholog_length,alignment_length,target_ortholog_identity,gap_openings,mismatches,e_value,bit_score,target_ortholog_coverage"
    diamond blastp \
        --query "${queries}" \
        --db "${db}" \
        --ultra-sensitive \
        --evalue 1e-3 \
        --outfmt 6 qseqid sseqid qlen slen length pident gapopen mismatch evalue bitscore \
        --max-target-seqs 25 \
        --threads ${task.cpus} \
    | sort -k6 -n \
    | awk -v OFS='\\t' -v header="\${header//,/\$'\\t'}" '
        BEGIN { print header } 
        { 
            split(\$1, target_id, "|"); 
            split(\$2, ortho_id, "|");
            \$8=(\$8/100);
            print target_id[2], ortho_id[2], \$0, \$5/\$3 
        }
    ' \
    > hits.tsv

    # rough histogram
    awk -F'\\t' -v OFS='\\t' '
        BEGIN { print "ortholog_uniprot_id", "ortholog_count" }
        (NR == 1) { for ( i=0; i<=NF; i++ ) a[\$i]=i }
        (NR > 1) { c[\$a["ortholog_uniprot_id"]]++ }
        END { for ( p in c ) print p, c[p] }
    ' hits.tsv \
    > hit_count.tsv

    """

}