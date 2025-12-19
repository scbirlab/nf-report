process umaps_of_rbh_matrix {

    tag "${id}"
    label 'big_mem'

    publishDir( 
        "${params.outputs}/targets/plots", 
        mode: 'copy',
        saveAs: { "${id}.${it}" }
    )

    input:
    tuple val( id ), path( matrix ), path( rowdata ), path( coldata )

    output:
    tuple val( id ), path( 'umap-taxonomy.{pdf,csv}' ), emit: strains
    tuple val( id ), path( 'umap-targets.{pdf,csv}' ), emit: targets
    tuple val( id ), path( 'umap-targets-taxon.{png,csv}' ), emit: targets_taxon


    script:
    """
    #!/usr/bin/env python
    from carabiner import print_err
    from carabiner.mpl import add_legend, grid, figsaver
    import pandas as pd
    import numpy as np
    from umap import UMAP

    figsave = figsaver(format="pdf", output_dir=".")

    rbh_m = pd.read_csv("${matrix}", sep="\\t", index_col=0)
    rbh_row_data = pd.read_csv("${rowdata}", sep="\\t", index_col=0)
    rbh_col_data = pd.read_csv("${coldata}", sep="\\t", index_col=0)
    rbh_m = pd.concat([rbh_m, rbh_row_data], axis=1)
    rbh_m = rbh_m.set_index(rbh_row_data.columns.tolist(), append=True)
    rbh_m.columns = pd.MultiIndex.from_frame(rbh_col_data.reset_index())
    print_err(rbh_m.head())

    reducer = UMAP(
        # n_neighbors=100,
        # min_dist=.5,
        random_state=42,
        # metric="cosine",
    )
    bacteria_embedding = pd.DataFrame(
        reducer.fit_transform(rbh_m.T),
        index=rbh_m.T.index,
        columns=["UMAP 1", "UMAP 2"],
    )

    target_weights_human = (
        rbh_m
        .groupby("target_is_human")
        .apply(
            lambda x: x.mean(axis=0),
        )
        .T
    )

    target_weights_bacteria = (
        rbh_m
        .groupby("target_is_bacteria")
        .apply(
            lambda x: x.mean(axis=0),
        )
        .T
    )

    fig, axes = grid(
        ncol=4, 
        aspect_ratio=1.2, 
        panel_size=4.,
    )
    axes[0].scatter(
        *bacteria_embedding.values.T,
        s=.1,
        color="lightgrey",
    )

    for tax_order, tax_data in bacteria_embedding.groupby("ortholog_order"):
        axes[1].scatter(
            tax_data.values[:,0],
            tax_data.values[:,1],
            s=.1,
            label=tax_order,
        )

    for ax, (w_name, w) in zip(
        axes[2:], 
        [
            ("human", target_weights_human),
            ("bacteria", target_weights_bacteria),
        ]
    ):
        sc = ax.scatter(
            *bacteria_embedding.values.T,
            s=.1,
            c=target_weights_human.iloc[:,1].values,
            cmap="magma",
            vmin=0., #vmax=1.,
        )
        fig.colorbar(sc, ax=ax, label=w_name)

    for ax in axes:
        ax.set(
            xlabel="UMAP 1",
            ylabel="UMAP 2",
        )
    figsave(fig, "umap-taxonomy", df=bacteria_embedding.reset_index())

    # ============================

    reducer = UMAP(
        # n_neighbors=50,
        min_dist=.5,
        # metric="cosine",
        random_state=42,
    )
    target_embedding = pd.DataFrame(
        reducer.fit_transform(rbh_m),
        index=rbh_m.index,
        columns=["UMAP 1", "UMAP 2"],
    )

    fig, axes = grid(ncol=9, aspect_ratio=1.2, panel_size=4.)
    axes[0].scatter(
        *target_embedding.values.T,
        s=.1,
        color="lightgrey"
    )

    for go_process, go_data in target_embedding.groupby("target_ec_number"):
        axes[1].scatter(
            go_data.values[:,0],
            go_data.values[:,1],
            s=.1,
            label=go_process,
        )
    axes[1].set(title="EC number")

    for tax_l1, tax_data in target_embedding.groupby("target_taxon_l1"):
        axes[2].scatter(
            tax_data.values[:,0],
            tax_data.values[:,1],
            s=.1,
            label=tax_l1,
        )
    axes[2].set(title="Taxon L1")
    add_legend(axes[2])

    for col, ax in zip(
        ("entropy", "sparsity", "mean_conservation", "median_conservation"), 
        axes[3:],
    ):
        sc = ax.scatter(
            *target_embedding.values.T,
            s=.1,
            c=rbh_row_data[col].loc[target_embedding.index.get_level_values("target_uniprot_id")].values,
            cmap="magma",
            vmin=0., vmax=1.,
        )
        fig.colorbar(sc, ax=ax, label=col)
        ax.set(title=col)


    for col, ax in zip(("target_is_human", "target_is_bacteria"), axes[7:]):
        ax.scatter(
            *target_embedding.values.T,
            s=.1,
            c="lightgrey",
        )
        this_m = target_embedding.query(col)
        ax.scatter(
            *this_m.values.T,
            s=.1,
            c="C1",
        )
        ax.set(title=col)
        
    for ax in axes:
        ax.set(
            xlabel="UMAP 1",
            ylabel="UMAP 2",
        )
    figsave(fig, "umap-targets", df=target_embedding.reset_index())

    # ============================
    from itertools import groupby

    taxonomy_weights = (
        rbh_m
        .query("target_taxon_l1 != 'Bacteria'").T
        .groupby(["ortholog_phylum", "ortholog_order"])
        .apply(
            lambda x: x.mean(axis=0),
        )
        .T
    )
    print_err(taxonomy_weights.shape)
    n_phyla = taxonomy_weights.columns.get_level_values("ortholog_phylum").nunique()
    grouped = groupby(
        sorted(
            taxonomy_weights.columns, 
            key=lambda x: x[0],
        ), 
        key=lambda x: x[0],
    )
    grouped = [
        (phylum, tuple(g)) 
        for phylum, g in grouped
    ]
    n_cols = max(len(g) for _, g in grouped)

    NCOL = 15
    fig, axes = grid(
        ncol=n_cols,
        nrow=n_phyla,
        aspect_ratio=1.1,
        panel_size=4.,
        squeeze=False,
    )

    for axrow, (phylum, taxa) in zip(axes, grouped):
        for ax, tax in zip(axrow, taxa):
            sc = ax.scatter(
                *target_embedding.query("target_taxon_l1 != 'Bacteria'").values.T,
                s=.1,
                color="lightgrey",
            )
            these_weights = taxonomy_weights[taxonomy_weights[tax] > .1][tax]
            this_tax_m = target_embedding.loc[these_weights.index.get_level_values("target_uniprot_id")]
            # print(this_tax_m.shape)
            sc = ax.scatter(
                *this_tax_m.values.T,
                s=.1,
                c=these_weights.values,
                cmap="magma",
                vmin=0., vmax=1.,
            )
            fig.colorbar(sc, ax=ax, label=tax)
            ax.set(title=" : ".join(tax))
        for ax in axrow[len(taxa):]:
            ax.set_axis_off()

    figsaver(format="png", output_dir=".")(fig, "umap-targets-taxon", df=target_embedding.reset_index())

    
    """
}
