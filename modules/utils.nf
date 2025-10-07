process stack_tables {

    tag "${id} -> ${directory}/${id}.${filename}"

    publishDir( 
        "${params.outputs}/${directory}", 
        mode: 'copy',
        saveAs: { (directory && filename) ? "${id}.${filename}" : null },
    )


    input:
    tuple val( id ), path( tables, stageAs: 'table-???.txt' )
    val directory
    val filename

    output:
    tuple val( id ), path( 'stacked.tsv' )

    script:
    """
    #!/usr/bin/env python

    from functools import partial
    from glob import glob
    import pandas as pd

    (
        pd.concat(
            map(
                partial(pd.read_csv, sep="\\t"),
                glob("table-???.txt"),
            ),
            axis=0,
        )
        .to_csv("stacked.tsv", sep="\\t", index=False)
    )
    
    """
}


process merge_tables {

    tag "${id}:${table1}:${table2}|${how} -> ${directory}/${id}.${filename}"

    errorStrategy 'retry'  // sometimes container fails to load
    maxRetries 2

    publishDir( 
        "${params.outputs}/${directory}", 
        mode: 'copy',
        saveAs: { (directory && filename) ? "${id}.${filename}" : null },
        // enabled: { directory && filename },
    )

    input:
    tuple val( id ), path( table1, stageAs: 'left.tsv' ), path( table2, stageAs: 'right.tsv' )
    val how
    val directory
    val filename

    output:
    tuple val( id ), path( 'merged.tsv' )

    script:
    """
    #!/usr/bin/env python

    import pandas as pd
    (
        pd.merge(
            pd.read_csv("${table1}", sep="\\t"),
            pd.read_csv("${table2}", sep="\\t"),
            how="${how}",
        )
        .drop_duplicates()
        .to_csv("merged.tsv", sep="\\t", index=False)
    )
    

    """
}


process split_csv {

    tag "${id}:${table}"

    // errorStrategy 'retry'  // sometimes container fails to load
    // maxRetries 2

    // publishDir( 
    //     "${params.outputs}/${directory}", 
    //     mode: 'copy',
    //     saveAs: { (directory && filename) ? "${id}.${filename}" : null },
    //     // enabled: { directory && filename },
    // )

    input:
    tuple val( id ), path( table )
    val chunksize
    // val how
    // val directory
    // val filename

    output:
    tuple val( id ), path( 'chunk-*.tsv' )

    script:
    """
    #!/usr/bin/env python

    import pandas as pd
    for i, chunk in enumerate(pd.read_csv("${table1}", sep="\\t", chunksize=${chunksize})):
        chunk.to_csv(f"chunk-{i}.tsv", sep="\\t", index=False)   

    """
}


process merge_tox_gnomad {

    tag "${table1}:${table2}:${table3}"

    publishDir( 
        "${params.outputs}/toxicity", 
        mode: 'copy',
    )

    input:
    path( table1 )
    path( table2 )
    path( table3 )
    val how

    output:
    path( 'target-tox-gnomad.tsv' )

    script:
    """
    #!/usr/bin/env python

    (
        pd.merge(
            (
                pd.read_csv("${table1}", sep="\\t")
                .query("target_taxon_id == 9606")
            ),
            pd.read_csv("${table2}", sep="\\t"),
            how="${how}",
        )
        .drop_duplicates()
        .merge(
            pd.read_csv("${table3}", sep="\\t"),
            how="${how}",
        )
        .drop_duplicates()
        .to_csv(
            "target-tox-gnomad.tsv", 
            sep="\\t", 
            index=False,
        )
    )
    
    """
}


process concat_files {

    tag "${id}"

    input:
    tuple val( id ), path( '*.txt' )

    output:
    tuple val( id ), path( 'concat.txt' )

    script:
    """
    cat *.txt > concat.txt
    
    """
}


process filter_target_list {

    tag "${id}: LOEUF ≤ ${min_loeuf}, ID > ${min_identity}, cov. > ${min_coverage}"

    publishDir( 
        "${params.outputs}/targets", 
        mode: 'copy',
        saveAs: { "${id}.${it}" }
    )

    input:
    tuple val( id ), path( table )
    val min_loeuf
    val min_identity
    val min_coverage

    output:
    tuple val( id ), path( 'conserved_hits.tsv' )

    script:
    """
    #!/usr/bin/env python
    import pandas as pd

    mammalia = "Mammalia"
    df = pd.read_csv("${table}", sep="\\t")
    print(df.head())

    potentially_toxic_targets = set(
        df
        .query(
            "target_taxon_id == 9606 "
            "and (LOEUF <= ${min_loeuf} or LOEUF.isna()) "
            "and target_ortholog_identity > ${min_identity} "
            "and target_ortholog_coverage > ${min_coverage}"
        )
        ["ortholog_uniprot_id"]
        .tolist()
    )

    (
        df
        .assign(
            potentially_toxic=lambda x: (
                x["ortholog_uniprot_id"]
                .isin(potentially_toxic_targets)
            )
        )
        .query(
            "(target_taxon_id == 9606 or target_taxon_l2 != @mammalia) "
            "and target_ortholog_identity > ${min_identity} "
            "and target_ortholog_coverage > ${min_coverage} "
        )
        .sort_values("target_ortholog_identity")
        .groupby(["target_uniprot_id", "target_name"])
        .tail(1)
        .groupby(["ortholog_uniprot_id", "target_taxon_id"])
        .tail(1)
        .drop_duplicates()
        .to_csv("conserved_hits.tsv", sep="\\t", index=False)
    )
    
    """
}