process factorise_nmf {

    tag "${id}"
    label 'big_mem'
    time '48h'

    publishDir( 
        "${params.outputs}/targets/nmf", 
        mode: 'copy',
        saveAs: { "${id}.${it}" }
    )

    input:
    tuple val( id ), path( matrix )

    output:
    tuple val( id ), path( '*.tsv' ), emit: factors
    tuple val( id ), path( 'plots/*.{csv,png}' ), emit: plots

    script:
    """
    #!/usr/bin/env python

    from collections import defaultdict
    import os
    import sys
    
    os.makedirs("mpl-tmp")
    os.environ["MPLCONFIGDIR"] = "mpl-tmp"
    
    from carabiner import print_err
    from carabiner.mpl import grid, figsaver
    import numpy as np
    import pandas as pd
    from sklearn.decomposition import NMF
    from tqdm.auto import tqdm

    figsave = figsaver(
        output_dir="plots",
        format="png",
        dpi=600,
    )
    

    def find_knee(x, y):
        # x, y already sorted by x
        x0, y0   = x[0],   y[0]          # first point
        x1, y1   = x[-1],  y[-1]         # last  point
        
        # vector form of the chord
        dx, dy   = x1 - x0, y1 - y0
        norm     = np.hypot(dx, dy)
        
        # perpendicular distance of every point to the chord
        dist = np.abs(dy*(x - x0) - dx*(y - y0)) / norm
        knee_idx = dist.argmax()
        return x[knee_idx], y[knee_idx]


    def nmf(X, ncomps=2, max_iter=1_000, seed=42, **kwargs):
        factorizer = NMF(
            n_components=int(ncomps), #'auto',
            max_iter=1_000,
            random_state=seed,
            **kwargs,
        )
        
        factorized = factorizer.fit(X)
        H, W = factorized.components_, factorizer.transform(X)
        m_approx = W @ H
        return factorizer, (W, H), factorized.reconstruction_err_


    def nmf_scan(X, num=10, max_iter=1_000, seed=42, **kwargs):

        n_comps = list(
            np.ceil(np.geomspace(2, min(max(1, *X.shape), 500) // 2, num=num))
            .astype(int)
        )
        print_err(f"Testing {n_comps=}")
        
        results = defaultdict(list)
        for n_comp in tqdm(n_comps):
            _, _, err = nmf(X, ncomps=n_comp)
            results["n_comp"].append(n_comp)
            results["err"].append(err)
        
        results = pd.DataFrame(results)
        knee_x, knee_y = find_knee(*results[["n_comp", "err"]].T.values)
        print_err(knee_x, knee_y)
        factorizer, (W, H), _ = nmf(X, ncomps=int(knee_x))
        return results, (knee_x, knee_y), factorizer, (W, H)

    m = pd.read_csv("${matrix}", sep="\\t", index_col=0)

    results, (knee_x, knee_y), factorizer, (W, H) = nmf_scan(
        m,
        num=100,    
    )
    m_approx = factorizer.inverse_transform(W)

    H = pd.DataFrame(
        H, 
        index=pd.Index(list(range(H.shape[0])), name="nmf_factors"),
        columns=pd.Index(m.columns.values, name="ortholog_taxon_id"),
    )
    H.to_csv(f"{H.columns.names[0]}.tsv", sep="\\t")
    W = pd.DataFrame(
        W, 
        index=pd.Index(m.index.values, name="target_uniprot_id"),
        columns=pd.Index(list(range(W.shape[1])), name="nmf_factors"),
    )
    W.to_csv(f"{W.index.names[0]}.tsv", sep="\\t")

    fig, axes = grid()
    axes.plot(
        "n_comp", "err",
        data=results,
    )
    axes.scatter(knee_x, knee_y, s=50.)
    axes.set(
        xlabel="Number of components",
        ylabel="Reconstruction error",
    )
    figsave(
        fig,
        "recon-err-trace",
        df=results,
    )

    modelled_vs_observed = pd.DataFrame(
        {
            "modelled": m_approx.ravel(),
            "observed": m.values.ravel(),
        },
    )
    modelled_vs_observed_sampled = modelled_vs_observed.sample(
        min(1_000_000, modelled_vs_observed.shape[0]), 
        random_state=42,
    )
    fig, axes = grid()
    axes.scatter(
        "modelled",
        "observed",
        data=modelled_vs_observed_sampled,
        s=1.,
    )
    axes.plot(
        axes.get_ylim(),
        axes.get_ylim(),
        color="lightgrey",
    )
    axes.set(
        title=f"Recon. error @ k = {H.shape[0]}",
        xlabel="Modelled",
        ylabel="Observed",
    )

    figsave(
        fig,
        "recon-err-scatter",
        df=modelled_vs_observed_sampled,
    )
    
    """
}
