process fetch_gnomad_constraints {

    tag "v${version}"

    publishDir( 
        "${params.outputs}/gnomad", 
        mode: 'copy',
        saveAs: { "${it}" }
    )

    input:
    val version

    output:
    path "*.constraint_metrics_.tsv"

    script:
    """
    curl https://storage.googleapis.com/gcp-public-data--gnomad/release/${version}/constraint/gnomad.v${version}.constraint_metrics.tsv > gnomad.v${version}.constraint_metrics.tsv

    python -c '
    import pandas as pd
    transcript = "protein_coding"
    (   
        pd.read_csv(
            "gnomad.v${version}.constraint_metrics.tsv", 
            sep="\\t",
        )
        .query("gene.str.len() > 0")
        .query("canonical and mane_select and transcript_type == @transcript")
        .assign(taxon_id=9606)
        [["taxon_id", "gene", "lof_hc_lc.pLI", "lof.oe_ci.upper"]]
        .rename(columns={
            "taxon_id": "target_taxon_id",
            "gene": "target_gene_symbol",
            "lof_hc_lc.pLI": "pLI", 
            "lof.oe_ci.upper": "LOEUF",
        })
        .sort_values("LOEUF")
        .groupby(["target_taxon_id", "target_gene_symbol"])
        .head(1)
        .to_csv(
            "gnomad.v${version}.constraint_metrics_.tsv", 
            sep="\\t", 
            index=False,
        )
    )
    '

    """

}