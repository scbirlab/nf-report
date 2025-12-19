process describe {

    tag "${id}"
    label 'big_mem'

    publishDir( 
        "${params.outputs}/orthology", 
        mode: 'copy',
        saveAs: { "${id}.${it}" }
    )

    input:
    tuple val( id ), path( table )

    output:
    tuple val( id ), path( 'plots/alignment-metrics.{csv,png}' ), emit: plots

    script:
    """
    #!/usr/bin/env python

    import os
    import sys
    
    os.makedirs("mpl-tmp")
    os.environ["MPLCONFIGDIR"] = "mpl-tmp"

    from carabiner.mpl import figsaver, scattergrid
    import pandas as pd
    import numpy as np


    figsave = figsaver(
        output_dir="plots",
        format="png",
        dpi=600,
    )

    HIST_VALUES = (
        "alignment_length",
        "mismatches",
        "gap_openings",
        "target_ortholog_identity",
        "target_ortholog_coverage",
        "e_value",
        "bit_score",
    )
    LOG_VALUES = (
        "e_value",
        "bit_score",
        "alignment_length",
        "mismatches",
        "gap_openings",
    )

    df = pd.read_csv("${table}", sep="\\t")
    print(f"{df['ortholog_taxon_id'].nunique()=}", file=sys.stderr)
    print(f"{df.shape=}", file=sys.stderr)

    df_sampled = df.sample(min(1_000_000, df.shape[0]), random_state=42)
    df_sampled = df_sampled.assign(**{col: df_sampled[col].astype(np.float64) for col in LOG_VALUES})
    fig, axes = scattergrid(
        df_sampled,
        grid_columns=list(HIST_VALUES),
        log=LOG_VALUES,
        n_bins=80,
        scatter_opts={"s": .1},
        hist_opts={"density": True},
    )
    figsave(fig, "alignment-metrics", df=df_sampled)


    """
}

process find_coverage_cutoff {

    tag "${id}"
    label 'big_mem'

    publishDir( 
        "${params.outputs}/orthology", 
        mode: 'copy',
        saveAs: { "${id}.${it}" }
    )

    input:
    tuple val( id ), path( full_table ), path( rbh_table )

    output:
    tuple val( id ), path( 'coverage-cutoff.txt' ), emit: cutoff
    tuple val( id ), path( 'plots/cov-cutoff*.{csv,png}' ), emit: plots

    script:
    """
    #!/usr/bin/env python

    import os
    import sys
    
    os.makedirs("mpl-tmp")
    os.environ["MPLCONFIGDIR"] = "mpl-tmp"
    
    from carabiner import print_err
    from carabiner.mpl import figsaver, grid, scattergrid
    import pandas as pd
    import numpy as np
    from sklearn.mixture import GaussianMixture
    from scipy.signal import argrelextrema

    figsave = figsaver(
        output_dir="plots",
        format="png",
        dpi=600,
    )

    HIST_VALUES = (
        "alignment_length",
        "mismatches",
        "gap_openings",
        "target_ortholog_identity",
        "target_ortholog_coverage",
        "e_value",
        "bit_score",
    )
    LOG_VALUES = (
        "e_value",
        "bit_score",
        "alignment_length",
        "mismatches",
        "gap_openings",
    )

    scattergrid_kwargs = {
        "grid_columns": HIST_VALUES,
        "grid_rows": ["target_ortholog_coverage", "target_ortholog_identity"],
        "log": LOG_VALUES,
        "n_bins": 80,
        "scatter_opts": {"s": .1},
        "hist_opts": {"density": True},
    }


    def find_coverage_nadir(coverage_values, resolution=10_000):
        xs = np.linspace(coverage_values.min(), coverage_values.max(), num=resolution)
        gmm = GaussianMixture(n_components=4).fit(coverage_values.reshape(-1, 1))
        logprob = gmm.score_samples(xs.reshape(-1,1))
        ys = np.exp(logprob)
        minima = argrelextrema(ys, np.less)[0]
        return (
            pd.DataFrame({"cutoff": xs, "logprob": logprob}), 
            (xs[minima[-1]], logprob[minima[-1]]),
        )

    df = pd.read_csv("${full_table}", sep="\\t")
    print(f"{df['ortholog_taxon_id'].nunique()=}", file=sys.stderr)
    print(f"{df.shape=}", file=sys.stderr)

    logprob, (COVERAGE_CUTOFF, cutoff_y) = find_coverage_nadir(df["target_ortholog_coverage"].values)
    print_err(f"{COVERAGE_CUTOFF=}")
    with open("coverage-cutoff.txt", "w") as f:
        print(COVERAGE_CUTOFF, file=f)

    fig, axes = grid()
    axes.plot(
        "cutoff",
        "logprob",
        data=logprob,
    )
    axes.scatter(
        [COVERAGE_CUTOFF],
        [cutoff_y],
        s=10.,
        c="C1",
    )
    axes.set(
        xlabel="Coverage cutoff",
        ylabel="Log prob. of minimum",
    )
    figsave(fig, "cov-cutoff-loglik", df=logprob)

    df_sampled = df.sample(min(1_000_000, df.shape[0]), random_state=42)
    df_sampled = df_sampled.assign(**{col: df_sampled[col].astype(np.float64) for col in LOG_VALUES})
    fig, axes = scattergrid(
        df_sampled,
        **scattergrid_kwargs,
    )
    figsave(fig, "cov-cutoff", df=df_sampled)
    for _x, ax in zip(HIST_VALUES, axes[0]):
        if _x == "target_ortholog_coverage":
            plotf = ax.axvline
            ax.plot(logprob["cutoff"], np.exp(logprob["logprob"]), color="dimgrey")
        else:
            plotf = ax.axhline
        plotf(COVERAGE_CUTOFF, color="lightgrey")
    figsave(fig, "cov-cutoff2", df=df_sampled)


    df_rbh = pd.read_csv("${rbh_table}", sep="\\t")
    rbh_sampled = df_rbh.sample(min(1_000_000, df_rbh.shape[0]))
    fig, axes = scattergrid(
        rbh_sampled,
        **scattergrid_kwargs,
    )

    for _x, ax in zip(HIST_VALUES, axes[0]):
        plotf = ax.axhline
        plotf(COVERAGE_CUTOFF, color="lightgrey")
    figsave(fig, "cov-cutoff-rbh", df=rbh_sampled)

    """

}
