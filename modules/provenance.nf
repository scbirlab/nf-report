process describe_provenance {

    tag "${id}"
    label 'big_mem'

    publishDir( 
        "${params.outputs}/orthology", 
        mode: 'copy',
        saveAs: { "${id}.${it}" }
    )

    input:
    tuple val( id ), path( table ), path( cutoff ), path( taxon_table )

    output:
    tuple val( id ), path( 'plots/{ortholog-yield-per-target,target-proportions,target-proportions-no-bacteria}.{csv,png}' ), emit: plots

    script:
    """
    #!/usr/bin/env python

    import os
    import sys
    
    os.makedirs("mpl-tmp")
    os.environ["MPLCONFIGDIR"] = "mpl-tmp"
    
    from carabiner.mpl import figsaver, grid, scattergrid
    import pandas as pd
    pd.options.mode.copy_on_write = True
    import numpy as np
    import seaborn as sns

    figsave = figsaver(
        output_dir="plots",
        format=["png", "svg", "pdf"],
        dpi=600,
    )

    def make_target_ortholog_summary(
        df, 
        query=None
    ):
        if query is not None and isinstance(query, str):
            df = df.query(query)
        target_ortholog_df = (
            df
            [[col for col in df if col.startswith(("target_", "ortholog_"))]]
            .drop_duplicates()
            .groupby(["target_taxon_id", "target_taxon_l1", "target_taxon_l2", "target_organism_name", "target_is_human", "target_is_bacteria"])
            [["target_uniprot_id", "target_ec_number", "ortholog_taxon_id", "ortholog_class", "ortholog_order", "ortholog_uniprot_id"]]
            .nunique()
        )
        return target_ortholog_df

    tax_df = (
        pd.read_csv("${taxon_table}", sep=",")
        .rename(columns={"organism_id": "ortholog_taxon_id"})
        .assign(ortholog_taxon_id=lambda x: x["ortholog_taxon_id"].astype(str))
    )
    tax_df = tax_df.rename(columns={
        col: f"ortholog_{col}" for col in tax_df 
        if not col.startswith("ortholog_")
    })
    df = pd.read_csv("${table}", sep="\\t").assign(ortholog_taxon_id=lambda x: x["ortholog_taxon_id"].astype(str))
    print(f"{df['ortholog_taxon_id'].nunique()=}", file=sys.stderr)
    print(f"{df.shape=}", file=sys.stderr)
    #if df.shape[0] > 50_000_000:
    #    df = df.sample(n=50_000_000, random_state=42)
    df = df.merge(tax_df, how="left")
    print(f"{df['ortholog_taxon_id'].nunique()=}", file=sys.stderr)
    print(f"{df.shape=}", file=sys.stderr)

    with open("${cutoff}", "r") as f:
        cutoff = f.readlines()[0].strip()
    cutoff = float(cutoff)
    print(f"{cutoff=}", file=sys.stderr)
    df["pass_coverage"] = df["target_ortholog_coverage"] >= cutoff
    print(f"{df.shape=}", file=sys.stderr)

    summary_all = make_target_ortholog_summary(df)
    summary_all["group"] = "All"
    print(f"{summary_all.shape=}", file=sys.stderr)
    summary_filtered = make_target_ortholog_summary(df, "pass_coverage")
    summary_filtered["group"] = "Cov. filter"
    print(f"{summary_filtered.shape=}", file=sys.stderr)
    target_ortholog_df = pd.concat([
        summary_all,
        summary_filtered,
    ], axis=0)
    print(f"{target_ortholog_df.shape=}", file=sys.stderr)

    fig, axes = scattergrid(
        target_ortholog_df,
        grid_columns=[col for col in target_ortholog_df if col.startswith("target_")],
        grid_rows=[col for col in target_ortholog_df if col.startswith("ortholog_")],
        log=target_ortholog_df.columns.tolist(),
        grouping="group",
        aspect_ratio=1.4,
    )
    figsave(fig, "ortholog-yield-per-target", df=target_ortholog_df)


    taxonomy_levels = (
        "target_taxon_l1", "target_taxon_l2", "target_is_human", "target_is_bacteria"
    )
    cats_to_count = (
        "target_uniprot_id",
        "ortholog_taxon_id",
        "ortholog_uniprot_id",
    )

    fig, axes = grid(
        ncol=target_ortholog_df["group"].nunique(), 
        nrow=len(taxonomy_levels), 
        aspect_ratio=2.5,
    )
    results = []
    for axcol, (gname, g_df) in zip(axes.T, target_ortholog_df.groupby("group")):
        for ax, tax in zip(axcol, taxonomy_levels):
            df_grouped = (
                g_df
                .groupby(tax)
                [list(cats_to_count)]
                .sum()
            )
            df_grouped_norm = df_grouped / df_grouped.sum(axis=0)
            df_grouped_norm = (
                df_grouped_norm
                .reset_index()
                .melt(id_vars=tax)
                .set_index(["variable", tax])
            )
            results.append(
                df_grouped
                .rename(columns={tax: "category_level"})
                .assign(category=tax)
            )
            df_grouped_norm.unstack().plot.barh(stacked=True, ax=ax)
            ax.set(ylabel="", title=f"{gname} : {tax}")
            sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    results = pd.concat(results, axis=0)
    figsave(fig, "target-proportions", df=results)

    
    taxonomy_levels = taxonomy_levels[:-1]
    fig, axes = grid(
        ncol=target_ortholog_df["group"].nunique(), 
        nrow=len(taxonomy_levels), 
        aspect_ratio=2.5,
    )
    results = []
    for axcol, (gname, g_df) in zip(
        axes.T, 
        target_ortholog_df.query("not target_is_bacteria").groupby("group"),
    ):
        for ax, tax in zip(axcol, taxonomy_levels):
            df_grouped = (
                g_df
                .groupby(tax)
                [list(cats_to_count)]
                .sum()
            )
            results.append(
                df_grouped
                .rename(columns={tax: "category_level"})
                .assign(category=tax)
            )
            df_grouped_norm = df_grouped / df_grouped.sum(axis=0)
            df_grouped_norm = (
                df_grouped_norm
                .reset_index()
                .melt(id_vars=tax)
                .set_index(["variable", tax])
            )
            df_grouped_norm.unstack().plot.barh(stacked=True, ax=ax)
            ax.set(ylabel="", title=f"{gname} : {tax}")
            sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    results = pd.concat(results, axis=0)
    figsave(fig, "target-proportions-no-bacteria", df=results)

    """
}
