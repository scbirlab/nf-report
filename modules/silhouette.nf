process silhouette {

    tag "${id}"
    label 'big_mem'
    time '3d'

    publishDir( 
        "${params.outputs}/targets/plots", 
        mode: 'copy',
        saveAs: { "${id}.${it}" }
    )

    input:
    tuple val( id ), path( matrix ), path( rowdata ), path( coldata )

    output:
    tuple val( id ), path( 'silhouette.{pdf,csv}' )


    script:
    """
    #!/usr/bin/env python
    from carabiner import print_err
    from carabiner.mpl import grid, figsaver
    import numpy as np
    import pandas as pd
    from sklearn.metrics import silhouette_score
    from tqdm.auto import tqdm

    BOOTSTRAP_N = 100
    figsave = figsaver(format="pdf", output_dir=".")
    rbh_m = pd.read_csv("${matrix}", sep="\\t", index_col=0)
    rbh_row_data = pd.read_csv("${rowdata}", sep="\\t", index_col=0)
    rbh_col_data = pd.read_csv("${coldata}", sep="\\t", index_col=0)
    rbh_m = pd.concat([rbh_m, rbh_row_data], axis=1)
    rbh_m = rbh_m.set_index(rbh_row_data.columns.tolist(), append=True)
    rbh_m.columns = pd.MultiIndex.from_frame(rbh_col_data.reset_index())
    print_err(rbh_m.head())

    rng = np.random.default_rng(seed=42)

    scores = []
    random_scores = {}
    p_vals = {}
    for level in tqdm(rbh_m.columns.names):
        print_err(level)
        this_m = rbh_m.T.query(f"not `{level}`.isna()")
        labels = this_m.index.get_level_values(level)
        print_err(labels)
        n_labels = len(set(labels))
        if n_labels > 2 and n_labels < (this_m.shape[0] - 2):
            this_silhouette_score = silhouette_score(
                this_m.values / this_m.values.sum(axis=1, keepdims=True),
                labels,
            )
            random_scores = np.asarray([
                silhouette_score(
                    rng.permuted(this_m, axis=0),
                    labels,
                ) for _ in range(BOOTSTRAP_N)
            ])
            scores.append({
                "level": level,
                "silhouette_score": this_silhouette_score,
                "bootstrap_median": np.median(random_scores),
                "bootstrap_low": np.median(random_scores) - np.percentile(random_scores, .05),
                "bootstrap_high": np.median(random_scores) - np.percentile(random_scores, .95),
                "p_val": max(1. / BOOTSTRAP_N, 1. - (this_silhouette_score > random_scores).mean()),
            })
    scores = pd.DataFrame(scores)
    print_err(scores)

    fig, axes = grid(panel_size=4.)
    axes.bar(
        "level",
        "silhouette_score",
        data=scores,
    )
    if scores.shape[0] > 0:
        axes.errorbar(
            "level",
            "bootstrap_median",
            np.stack([
                scores["bootstrap_low"].values, 
                scores["bootstrap_high"].values,
            ], axis=0),
            data=scores,
            color="lightgrey",
        )
    axes.set(
        xlabel="Level",
        ylabel="Silhouette score",
    )
    figsave(fig, "silhouette", df=scores)
    
    """
}
