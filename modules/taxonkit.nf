process fetch_taxonkit_db {

    tag "${url}"

    publishDir( 
        "${params.outputs}/taxonkit", 
        mode: 'copy',
    )

    input:
    val url

    output:
    path '.taxonkit'

    script:
    """
    set -euxo pipefail
    
    curl -v "${url}" -o taxdump.tar.gz
    tar -zxvf taxdump.tar.gz

    mkdir .taxonkit
    mv names.dmp nodes.dmp delnodes.dmp merged.dmp .taxonkit
    rm *.{dmp,prt,txt}

    """
}


process fetch_taxonomic_ranks {

    tag "${table}"

    publishDir( 
        "${params.outputs}/taxonkit", 
        mode: 'copy',
    )

    input:
    tuple path( table ), path( taxonkit_db )
    val column

    output:
    path "taxonomy.csv"

    script:
    """
    set -euox pipefail

    tax_id_col=\$(head -n1 "${table}" | tr , \$'\\n' | grep -nxF "${column}" | cut -d: -f1)

    head -n1 "${table}" \
    | awk -F, -v OFS=, -v col="\$tax_id_col" '
        { 
            print \$col, "domain", "kingdom", "phylum", "class", "order", "family", "genus", "species", "subspecies", "strain"
        }
        ' \
    > "taxonomy.csv"

    tail -n+2 "${table}" \
    | cut -f2 -d, \
    | sort -u -n \
    | taxonkit reformat \
        --data-dir "${taxonkit_db}" \
        --threads "${task.cpus}" \
        -I 1 \
        --pseudo-strain \
        -f '{d},{K},{p},{c},{o},{f},{g},{s},{S},{t}' \
    | tr \$'\\t' , \
    >> "taxonomy.csv"

    """
}