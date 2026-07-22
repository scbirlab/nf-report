process stack_tables {

    tag "${id} -> ${directory}/${id}.${filename}.gz"
    label 'big_mem'
    time '1d'
    stageInMode 'link'

    publishDir( 
        "${params.outputs}/${directory}", 
        mode: 'copy',
        saveAs: { (directory && filename) ? "${id}.${filename}.gz" : null },
    )


    input:
    tuple val( id ), path( tables, stageAs: 'table-?????/*' )
    val directory
    val filename

    output:
    tuple val( id ), path( 'stacked.tsv.gz' )

    script:
    """
    #!/usr/bin/env python

    from functools import partial
    from glob import glob
    import gzip

    import pandas as pd

    reader = partial(pd.read_csv, sep="\\t")
    files = glob("table-*/*")

    trial = {f: reader(f, nrows=2) for f in files}
    non_empty = [f for f, df in trial.items() if df.shape[0] > 0]
    if len(non_empty) > 0:
        heads = pd.concat([
            trial[f] for f in non_empty
        ], axis=0)
        all_cols = heads.columns

        with gzip.open("stacked.tsv.gz", "wb") as out:
            for i, f in enumerate(non_empty):
                df = reader(f).reindex(all_cols, axis=1)
                df.to_csv(
                    out, 
                    sep="\\t", 
                    index=False,
                    header=i == 0,
                )
    else:
        pd.DataFrame().to_csv("stacked.tsv.gz", sep="\\t", index=False)
    
    """
}


process subset_table {

    tag "${id}:${column} -> ${directory}/${id}.${filename}.gz"
    label 'big_mem'
    stageInMode 'link'

    publishDir( 
        "${params.outputs}/${directory}", 
        mode: 'copy',
        saveAs: { (directory && filename) ? "${id}.${filename}.gz" : null },
    )


    input:
    tuple val( id ), path( tables, stageAs: 'table-?????/*' )
    val column
    val directory
    val filename

    output:
    tuple val( id ), path( 'subset.tsv.gz' )

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
                glob("table-*/*"),
            ),
            axis=0,
        )
        .query("${column} == '${id}'")
        .to_csv("subset.tsv.gz", sep="\\t", index=False)
    )
    
    """
}


process merge_tables {

    tag "${id}|${how} -> ${directory}/${id}.${filename}.gz"
    stageInMode 'link'

    errorStrategy 'retry'  // sometimes container fails to load
    maxRetries 2

    publishDir( 
        "${params.outputs}/${directory}", 
        mode: 'copy',
        saveAs: { (directory && filename) ? "${id}.${filename}.gz" : null },
    )

    input:
    tuple val( id ), path( table1, stageAs: 'left-??/*' ), path( table2, stageAs: 'right-??/*' )
    val how
    val directory
    val filename

    output:
    tuple val( id ), path( 'merged.tsv.gz' )

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
        .to_csv("merged.tsv.gz", sep="\\t", index=False)
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

    output:
    tuple val( id ), path( 'chunk-*.tsv.gz' )

    script:
    """
    #!/usr/bin/env python

    import pandas as pd
    for i, chunk in enumerate(pd.read_csv("${table1}", sep="\\t", chunksize=${chunksize})):
        chunk.to_csv(f"chunk-{i}.tsv.gz", sep="\\t", index=False)   

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
    path( 'target-tox-gnomad.tsv.gz' )

    script:
    """
    #!/usr/bin/env python
    import pandas as pd

    (
        pd.read_csv("${table1}", sep="\\t")
        .query("target_taxon_id == 9606")
        .merge(
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
            "target-tox-gnomad.tsv.gz", 
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

    tag "${id}: LOEUF ≤ ${min_loeuf}, ID > ${min_identity}"

    publishDir( 
        "${params.outputs}/targets/conserved-hits", 
        mode: 'copy',
        saveAs: { "${id}.${it}" }
    )

    input:
    tuple val( id ), path( table ), path( min_coverage )
    val min_loeuf
    val min_identity

    output:
    tuple val( id ), path( 'conserved_hits.tsv.gz' )

    script:
    """
    #!/usr/bin/env python
    import pandas as pd

    with open("${min_coverage}", "r") as f:
        COVERAGE_CUTOFF = float(f.readlines()[0])

    mammalia = "Mammalia"
    df = pd.read_csv("${table}", sep="\\t")
    print(df.head())

    potentially_toxic_targets = set(
        df
        .query(
            "target_taxon_id == 9606 "
            "and (LOEUF <= ${min_loeuf} or LOEUF.isna()) "
            "and target_ortholog_identity > ${min_identity} "
            "and target_ortholog_coverage > @COVERAGE_CUTOFF"
        )
        ["ortholog_uniprot_id"]
        .unique()
    )

    (
        df
        .assign(
            potentially_toxic=lambda x: (
                x["ortholog_uniprot_id"]
                .isin(potentially_toxic_targets)
            ),
            ortholog_taxon_id="${id}",
        )
        .sort_values("bit_score")
        .groupby(["target_taxon_id", "target_uniprot_id", "ortholog_taxon_id"])
        .tail(1)
        .groupby(["target_taxon_id", "ortholog_taxon_id", "ortholog_uniprot_id"])
        .tail(1)
        .drop_duplicates()
        .query(
            "(target_taxon_id == 9606 or target_taxon_l2 != @mammalia) "
            "and target_ortholog_identity > ${min_identity} "
            "and target_ortholog_coverage > @COVERAGE_CUTOFF "
        )
        .to_csv("conserved_hits.tsv.gz", sep="\\t", index=False)
    )
    
    """
}

process make_rbh_matrix {

    tag "${id}"
    label 'big_mem'

    publishDir( 
        "${params.outputs}/targets", 
        mode: 'copy',
        saveAs: { "${id}.${it}" }
    )

    input:
    tuple val( id ), path( table ), path( taxonomy )

    output:
    tuple val( id ), path( 'rbh.tsv.gz' ), emit: table
    tuple val( id ), path( 'rbh_m.tsv.gz' ), emit: matrix
    tuple val( id ), path( 'rbh_m_coverage.tsv.gz' ), emit: matrix_coverage
    tuple val( id ), path( 'rbh_m.rowdata.tsv.gz' ), emit: row_data
    tuple val( id ), path( 'rbh_m.coldata.tsv.gz' ), emit: col_data

    script:
    """
    #!/usr/bin/env python
    import numpy as np
    import pandas as pd
    from scipy.stats import entropy


    def _entropy(x, axis=0, unnorm=False):
        axis_kwargs = {"axis": axis, "keepdims": True}
        x = np.asarray(x)
        x = x / x.sum(**axis_kwargs)
        H =  entropy(x, **axis_kwargs)
        return H / np.log2(x.shape[axis]) if not unnorm else H


    def sparsity_index(x, axis=0):
        axis_kwargs = {"axis": axis, "keepdims": True}
        x = np.asarray(x)
        x = x / x.sum(**axis_kwargs)
        n = x.shape[axis]
        sqrt_n = np.sqrt(n)
        return (
            sqrt_n - np.linalg.norm(x, 1, **axis_kwargs) 
            / np.linalg.norm(x, 2, **axis_kwargs)
        ) / (sqrt_n - 1.)
        


    df = pd.read_csv("${table}", sep="\\t")
    print(df.head())

    rbh = (
        df
        .sort_values("bit_score")
        .groupby(["target_taxon_id", "target_uniprot_id", "ortholog_taxon_id"])
        .tail(1)
        .groupby(["target_taxon_id", "ortholog_taxon_id", "ortholog_uniprot_id"])
        .tail(1)
    )
    rbh.to_csv("rbh.tsv.gz", sep="\\t", index=False)
    rbh_m = (
        rbh
        .assign(ortholog_taxon_id=lambda x: x["ortholog_taxon_id"].astype(str))
        .groupby(["target_uniprot_id", "ortholog_taxon_id"])
        .tail(1)
        .drop_duplicates()
        .pivot(
            index="target_uniprot_id",
            columns="ortholog_taxon_id",
            values="target_ortholog_identity",
        )
        .fillna(0.)
    )
    rbh_m.to_csv("rbh_m.tsv.gz", sep="\\t", index=True)
    rbh_m_cov = (
        rbh
        .assign(ortholog_taxon_id=lambda x: x["ortholog_taxon_id"].astype(str))
        .groupby(["target_uniprot_id", "ortholog_taxon_id"])
        .tail(1)
        .drop_duplicates()
        .pivot(
            index="target_uniprot_id",
            columns="ortholog_taxon_id",
            values="target_ortholog_coverage",
        )
        .fillna(0.)
    )
    rbh_m_cov.to_csv("rbh_m_coverage.tsv.gz", sep="\\t", index=True)

    tax_df = (
        pd.read_csv("${taxonomy}")
        .rename(columns={"organism_id": "taxon_id"})
    )
    rbh_col_data = (
        tax_df
        .rename(columns={col: f"ortholog_{col}" for col in tax_df})
        .assign(ortholog_taxon_id=lambda x: x["ortholog_taxon_id"].astype(str))
        .set_index("ortholog_taxon_id")
        .loc[rbh_m.columns.get_level_values("ortholog_taxon_id")]
    )
    rbh_col_data.to_csv("rbh_m.coldata.tsv.gz", sep="\\t")
    rbh_row_data = (
        df
        [[col for col in df if col.startswith("target_") and not "_ortholog_" in col]]
        .drop_duplicates()
        .groupby("target_uniprot_id")
        .tail(1)
        .set_index("target_uniprot_id")
        .loc[rbh_m.index]
        .assign(
            entropy=_entropy(rbh_m, axis=1),
            sparsity=sparsity_index(rbh_m, axis=1),
            mean_conservation=rbh_m.values.mean(axis=1),
            median_conservation=np.median(rbh_m.values, axis=1),
        )
    )
    rbh_row_data = (
        rbh_row_data
        .assign(
            entropy=_entropy(rbh_m, axis=1),
            sparsity=sparsity_index(rbh_m, axis=1),
            mean_conservation=rbh_m.values.mean(axis=1),
            median_conservation=np.median(rbh_m.values, axis=1),
        )
    )
    rbh_row_data.to_csv("rbh_m.rowdata.tsv.gz", sep="\\t")
    
    """
}
