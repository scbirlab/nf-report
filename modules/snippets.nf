/*
 * Python snippets shared by several processes. Each function returns code that
 * is interpolated into a process script (indented to match its callers).
 */

def _indented( String code, int spaces = 4 ) {
    // drop trailing whitespace first so the closing quotes' indent doesn't count
    return code.replaceAll( /\s+$/, "" ).stripIndent().trim().readLines().join( "\n" + ( " " * spaces ) )
}

// Columns shown in the alignment-metric plots, and the ones plotted on a log scale
def plot_constants_py() {
    return _indented("""
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
    """)
}

// Best reciprocal hits: best ortholog per target per taxon, then best target per ortholog
def rbh_py() {
    return _indented("""
        def best_reciprocal_hits(df):
            return (
                df
                .sort_values("bit_score")
                .groupby(["target_taxon_id", "target_uniprot_id", "ortholog_taxon_id"])
                .tail(1)
                .groupby(["target_taxon_id", "ortholog_taxon_id", "ortholog_uniprot_id"])
                .tail(1)
            )
    """)
}
