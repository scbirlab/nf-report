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
        format=["png", "svg", "pdf"],
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
    tuple val( id ), path( 'plots/rbh-before-after-hist.{csv,pdf}' ), emit: histograms

    script:
    """
    #!/usr/bin/env python

    import os
    import sys
    
    os.makedirs("mpl-tmp")
    os.environ["MPLCONFIGDIR"] = "mpl-tmp"
    
    from carabiner import print_err
    from carabiner.mpl import add_legend, figsaver, grid, scattergrid
    import pandas as pd
    import numpy as np
    from sklearn.mixture import GaussianMixture
    from scipy.signal import argrelextrema

    figsave = figsaver(
        output_dir="plots",
        format=["png", "svg", "pdf"],
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

    RBH_KEYS = [
        "target_taxon_id", 
        "target_uniprot_id", 
        "ortholog_taxon_id", 
        "ortholog_uniprot_id",
    ]

    scattergrid_kwargs = {
        "grid_columns": HIST_VALUES,
        "grid_rows": ["target_ortholog_coverage", "target_ortholog_identity"],
        "log": LOG_VALUES,
        "n_bins": 80,
        "scatter_opts": {"s": .1},
        "hist_opts": {"density": True},
    }


    COVERAGE_CUTOFF_FALLBACK = float("${params.coverage_cutoff_fallback}")


    def find_coverage_nadir(coverage_values, resolution=10_000):
        # Returns (density table, (cutoff, log-density at cutoff), method).
        # method is "inferred" only if a real interior minimum of the fitted density exists;
        # otherwise the configured fallback is used and labelled as such.
        coverage_values = np.asarray(coverage_values, dtype=float)
        coverage_values = coverage_values[np.isfinite(coverage_values)]
        if coverage_values.shape[0] < 4 or np.ptp(coverage_values) == 0:
            print_err("Too few or constant coverage values to fit a density; using fallback cutoff")
            return None, (COVERAGE_CUTOFF_FALLBACK, np.nan), "fallback_insufficient_data"
        xs = np.linspace(coverage_values.min(), coverage_values.max(), num=resolution)
        gmm = GaussianMixture(n_components=4, random_state=42).fit(coverage_values.reshape(-1, 1))
        logprob = gmm.score_samples(xs.reshape(-1, 1))
        # minima on the log density (monotone in the density, but no underflow)
        minima = argrelextrema(logprob, np.less)[0]
        density = pd.DataFrame({"cutoff": xs, "logprob": logprob})
        if len(minima) == 0:
            print_err("No local minimum in coverage density; using fallback cutoff")
            return density, (COVERAGE_CUTOFF_FALLBACK, np.nan), "fallback_no_minimum"
        return density, (xs[minima[-1]], logprob[minima[-1]]), "inferred"

    df = pd.read_csv("${full_table}", sep="\\t")
    print(f"{df['ortholog_taxon_id'].nunique()=}", file=sys.stderr)
    print(f"{df.shape=}", file=sys.stderr)

    logprob, (COVERAGE_CUTOFF, cutoff_y), CUTOFF_METHOD = find_coverage_nadir(df["target_ortholog_coverage"].values)
    print_err(f"{COVERAGE_CUTOFF=} ({CUTOFF_METHOD})")
    with open("coverage-cutoff.txt", "w") as f:
        print(COVERAGE_CUTOFF, file=f)
    os.makedirs("plots", exist_ok=True)
    pd.DataFrame({
        "cutoff": [COVERAGE_CUTOFF], 
        "method": [CUTOFF_METHOD],
    }).to_csv("plots/cov-cutoff-method.csv", index=False)

    ## Plot profile likelihood (only if a density could be fitted)
    if logprob is not None:
        fig, axes = grid()
        axes.plot(
            "cutoff",
            "logprob",
            data=logprob,
        )
        if np.isfinite(cutoff_y):
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

    ## Plot full grid of interactions
    df_sampled = df.sample(min(1_000_000, df.shape[0]), random_state=42)
    df_sampled = df_sampled.assign(**{col: df_sampled[col].astype(np.float64) for col in LOG_VALUES})
    fig, axes = scattergrid(
        df_sampled,
        **scattergrid_kwargs,
    )
    figsave(fig, "cov-cutoff", df=df_sampled)

    ## Plot smaller grid of interactions
    for _x, ax in zip(HIST_VALUES, axes[0]):
        if _x == "target_ortholog_coverage":
            plotf = ax.axvline
            if logprob is not None:
                ax.plot(logprob["cutoff"], np.exp(logprob["logprob"]), color="dimgrey")
        else:
            plotf = ax.axhline
        plotf(COVERAGE_CUTOFF, color="lightgrey")
    figsave(fig, "cov-cutoff2", df=df_sampled)


    ## Plot smaller grid of interactions after RBH
    df_rbh = pd.read_csv("${rbh_table}", sep="\\t")
    rbh_sampled = df_rbh.sample(min(1_000_000, df_rbh.shape[0]), random_state=42)
    fig, axes = scattergrid(
        rbh_sampled,
        **scattergrid_kwargs,
    )

    for _x, ax in zip(HIST_VALUES, axes[0]):
        plotf = ax.axvline if _x == "target_ortholog_coverage" else ax.axhline
        plotf(COVERAGE_CUTOFF, color="lightgrey")
    figsave(fig, "cov-cutoff-rbh", df=rbh_sampled)


    ## Plot historgram of identity and coverage before and after RBH
    fig, axes = grid(ncol=len(HIST_VALUES), aspect_ratio=1.2)

    for ax, col in zip(axes, HIST_VALUES):
        for label, _df in zip(["full", "rbh"], [df, df_rbh]):
            ax.hist(
                col,
                data=_df,
                bins=80,
                density=True,
                label=label,
                histtype="stepfilled",
                edgecolor="dimgrey",
                linewidth=1.,
                alpha=.7,
            )
        ax.set(
            xlabel=col,
            ylabel="Density",
        )
    add_legend(ax)
    figsaver(
        output_dir="plots",
        format="pdf",
        dpi=600,
    )(
        fig, 
        "rbh-before-after-hist", 
        df=df
            .merge(
                df_rbh[RBH_KEYS].drop_duplicates().assign(retained_in_rbh=True),
                on=RBH_KEYS,
                how="left",
            )
            .assign(retained_in_rbh=lambda x: x["retained_in_rbh"].fillna(False).astype(bool))
        ,
    )

    """

}
