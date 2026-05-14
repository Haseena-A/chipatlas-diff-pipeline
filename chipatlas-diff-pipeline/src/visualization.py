"""
Publication-quality visualisations for ChIP-Atlas differential analysis.

Generates four standard plots saved as both PNG (300 DPI) and SVG:

1. **Volcano plot**   – logFC vs −log₁₀(Q-value)
2. **Chromosome distribution** – significant region counts per chromosome
3. **Region size histogram** – log-scale size distribution by direction
4. **MA plot** – log₂ mean count vs logFC

Uses ``plotnine`` with the ``plotnine-prism`` theme for publication aesthetics.
Gracefully degrades to PNG-only if SVG export fails.
"""

from __future__ import annotations

import logging
import os
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# Colour palette (red = A, blue = B, grey = NS)
COLOUR_A    = "#E74C3C"
COLOUR_B    = "#3498DB"
COLOUR_NS   = "#95A5A6"


def generate_all_plots(
    results: Dict,
    output_dir: str = "diff_analysis_results",
) -> Optional[str]:
    """
    Generate all four standard visualisations and save to *output_dir*.

    Parameters
    ----------
    results : dict
        Output of :func:`~src.diff_analysis.run_diff_workflow`.
    output_dir : str
        Destination directory (created if absent).

    Returns
    -------
    str or None
        Path to the volcano plot PNG, or None if no data to plot.
    """
    try:
        from plotnine import (
            aes, element_text,
            geom_col, geom_histogram, geom_hline, geom_point, geom_vline,
            ggplot, labs, position_dodge,
            scale_color_manual, scale_fill_manual,
            scale_x_discrete, scale_x_log10, theme,
        )
        from plotnine_prism import theme_prism
    except ImportError as exc:
        logger.warning(
            "plotnine / plotnine-prism not installed – skipping plots (%s)", exc
        )
        print(f"   Warning: plotting skipped (install plotnine and plotnine-prism): {exc}")
        return None

    os.makedirs(output_dir, exist_ok=True)
    print("\n" + "=" * 70)
    print("GENERATING VISUALISATIONS")
    print("=" * 70 + "\n")

    df      = results["diff_regions"]
    params  = results["parameters"]
    desc_a  = params.get("description_a", "Group A")
    desc_b  = params.get("description_b", "Group B")
    genome  = params.get("genome", "hg38")
    qc_warn = results.get("qc_warnings", [])

    if len(df) == 0:
        print("   Warning: No differential regions to plot.")
        return None

    # 1. Volcano
    print("1. Generating volcano plot …")
    p1 = _plot_volcano(df, desc_a, desc_b)
    _save_plot(p1, output_dir, "volcano_plot")

    # 2. Chromosome distribution
    print("2. Generating chromosome distribution plot …")
    p2 = _plot_chromosome_distribution(df, desc_a, desc_b, qc_warn, genome)
    _save_plot(p2, output_dir, "chromosome_distribution", width=10, height=6)

    # 3. Region size
    print("3. Generating region size distribution …")
    p3 = _plot_region_size(df, desc_a, desc_b)
    _save_plot(p3, output_dir, "region_size_distribution")

    # 4. MA plot
    print("4. Generating MA plot …")
    p4 = _plot_ma(df, desc_a, desc_b)
    _save_plot(p4, output_dir, "ma_plot")

    print("\n✓ All visualisations generated successfully!")
    print("=" * 70 + "\n")
    return os.path.join(output_dir, "volcano_plot.png")


# ─────────────────────────── Plot builders ─────────────────────────────────── #


def _plot_volcano(df: pd.DataFrame, desc_a: str, desc_b: str):
    """Volcano plot: logFC vs −log₁₀(Q-value)."""
    from plotnine import (
        aes, geom_hline, geom_point, geom_vline,
        ggplot, labs, scale_color_manual, theme,
    )
    from plotnine_prism import theme_prism

    plot_df = df.copy()
    plot_df["neg_log10_q"] = -np.log10(plot_df["qvalue"].clip(lower=1e-300))
    sig = plot_df.get("significant", plot_df["qvalue"] < 0.05)
    plot_df["group"] = np.where(
        ~sig,
        "Not significant",
        np.where(
            plot_df["logFC"] > 0,
            f"{desc_a} enriched",
            f"{desc_b} enriched",
        ),
    )
    colours = {
        "Not significant":   COLOUR_NS,
        f"{desc_a} enriched": COLOUR_A,
        f"{desc_b} enriched": COLOUR_B,
    }
    return (
        ggplot(plot_df, aes(x="logFC", y="neg_log10_q", color="group"))
        + geom_point(size=1.5, alpha=0.6)
        + geom_hline(
            yintercept=-np.log10(0.05),
            linetype="dashed", color="gray", size=0.5,
        )
        + geom_vline(xintercept=0, linetype="solid", color="gray", size=0.3, alpha=0.5)
        + scale_color_manual(values=colours)
        + labs(
            title="Volcano Plot",
            x="Log₂ Fold Change",
            y="−log₁₀(Q-value)",
            color="Direction",
        )
        + theme_prism()
        + theme(legend_position="bottom")
    )


def _plot_chromosome_distribution(
    df: pd.DataFrame,
    desc_a: str,
    desc_b: str,
    qc_warnings: List[Dict],
    genome: str,
):
    """Grouped bar chart of significant regions per standard chromosome."""
    from plotnine import (
        aes, geom_col, ggplot, labs,
        position_dodge, scale_fill_manual, scale_x_discrete, theme,
    )
    from plotnine_prism import theme_prism
    from src.validation import get_standard_chroms

    sig_df = df[df.get("significant", df["qvalue"] < 0.05)] if len(df) else df
    standard = get_standard_chroms(genome)
    std_sig = sig_df[sig_df["chrom"].isin(standard)]
    n_other = len(sig_df) - len(std_sig)

    sex_confounds: set = set()
    for w in qc_warnings:
        if "sex chromosome confound" in w.get("issue", ""):
            c = w["issue"].split(" ")[0]
            sex_confounds.add(c.replace("chr", ""))

    chrom_order = _chrom_sort_key(std_sig["chrom"].unique())
    rows = []
    for c in chrom_order:
        label = c.replace("chr", "")
        if label in sex_confounds:
            label += "*"
        na = int(((std_sig["chrom"] == c) & (std_sig["direction"] == "A_enriched")).sum())
        nb = int(((std_sig["chrom"] == c) & (std_sig["direction"] == "B_enriched")).sum())
        if na or nb:
            rows.append({"chrom": label, "count": na, "group": f"{desc_a} enriched"})
            rows.append({"chrom": label, "count": nb, "group": f"{desc_b} enriched"})

    if not rows:
        from plotnine import ggplot, aes, labs
        return (
            ggplot(pd.DataFrame({"x": [0], "y": [0]}), aes("x", "y"))
            + labs(title="No significant differential regions by chromosome")
            + theme_prism()
        )

    plot_df = pd.DataFrame(rows)
    active  = [
        c.replace("chr", "") + ("*" if c.replace("chr", "") in sex_confounds else "")
        for c in chrom_order
        if c.replace("chr", "") in plot_df["chrom"].values
        or c.replace("chr", "") + "*" in plot_df["chrom"].values
    ]

    title = "Differential Regions by Chromosome (FDR < 0.05)"
    if n_other:
        title += f"\n({n_other} regions on non-standard contigs not shown)"
    if sex_confounds:
        title += "\n(* probable sex-chromosome artefact)"

    return (
        ggplot(plot_df, aes(x="chrom", y="count", fill="group"))
        + geom_col(position=position_dodge(preserve="single"), alpha=0.8)
        + scale_fill_manual(values={f"{desc_a} enriched": COLOUR_A, f"{desc_b} enriched": COLOUR_B})
        + scale_x_discrete(limits=active)
        + labs(title=title, x="Chromosome", y="Number of Significant Regions", fill="Direction")
        + theme_prism()
        + theme(legend_position="bottom")
    )


def _plot_region_size(df: pd.DataFrame, desc_a: str, desc_b: str):
    """Log-scale histogram of region sizes split by enrichment direction."""
    from plotnine import (
        aes, geom_histogram, ggplot, labs,
        scale_fill_manual, scale_x_log10, theme,
    )
    from plotnine_prism import theme_prism

    if "region_size" not in df.columns or df["region_size"].dropna().empty:
        return (
            ggplot(pd.DataFrame({"x": [0], "y": [0]}), aes("x", "y"))
            + labs(title="No region size data available")
            + theme_prism()
        )

    plot_df = df[["region_size", "direction"]].dropna().copy()
    plot_df = plot_df[plot_df["region_size"] > 0]
    plot_df["group"] = np.where(
        plot_df["direction"] == "A_enriched",
        f"{desc_a} enriched",
        f"{desc_b} enriched",
    )
    return (
        ggplot(plot_df, aes(x="region_size", fill="group"))
        + geom_histogram(bins=30, alpha=0.8, position="dodge")
        + scale_x_log10()
        + scale_fill_manual(values={f"{desc_a} enriched": COLOUR_A, f"{desc_b} enriched": COLOUR_B})
        + labs(
            title="Region Size Distribution",
            x="Region Size (bp, log scale)",
            y="Frequency",
            fill="Direction",
        )
        + theme_prism()
        + theme(legend_position="bottom")
    )


def _plot_ma(df: pd.DataFrame, desc_a: str, desc_b: str):
    """MA-style plot: log₂ mean normalised count vs logFC."""
    from plotnine import (
        aes, geom_hline, geom_point, ggplot, labs,
        scale_color_manual, theme,
    )
    from plotnine_prism import theme_prism

    if "mean_count_a" not in df.columns or "mean_count_b" not in df.columns:
        return (
            ggplot(pd.DataFrame({"x": [0], "y": [0]}), aes("x", "y"))
            + labs(title="No count data for MA plot")
            + theme_prism()
        )

    plot_df = df.copy()
    mean_expr = (plot_df["mean_count_a"] + plot_df["mean_count_b"]) / 2
    plot_df["log2_mean_count"] = np.log2(mean_expr.clip(lower=0.1))
    sig = plot_df.get("significant", plot_df["qvalue"] < 0.05)
    plot_df["group"] = np.where(
        ~sig,
        "Not significant",
        np.where(
            plot_df["logFC"] > 0,
            f"{desc_a} enriched",
            f"{desc_b} enriched",
        ),
    )
    colours = {
        "Not significant":    COLOUR_NS,
        f"{desc_a} enriched": COLOUR_A,
        f"{desc_b} enriched": COLOUR_B,
    }
    return (
        ggplot(plot_df, aes(x="log2_mean_count", y="logFC", color="group"))
        + geom_point(size=1.2, alpha=0.5)
        + geom_hline(yintercept=0, linetype="solid", color="gray", size=0.3, alpha=0.5)
        + scale_color_manual(values=colours)
        + labs(
            title="MA Plot",
            x="Log₂ Mean Normalised Count",
            y="Log₂ Fold Change",
            color="Direction",
        )
        + theme_prism()
        + theme(legend_position="bottom")
    )


# ─────────────────────────── Helpers ───────────────────────────────────────── #


def _chrom_sort_key(chromosomes) -> List[str]:
    """Sort chromosomes in natural genomic order."""
    def key(c: str):
        name = c.replace("chr", "")
        return (0, int(name)) if name.isdigit() else (1, name)
    return sorted(chromosomes, key=key)


def _save_plot(plot, output_dir: str, base_name: str, width: int = 8, height: int = 6) -> None:
    """Save a plotnine plot as PNG (300 DPI) and optionally SVG."""
    png_path = os.path.join(output_dir, f"{base_name}.png")
    plot.save(png_path, dpi=300, width=width, height=height, verbose=False)
    print(f"   Saved: {png_path}")

    svg_path = os.path.join(output_dir, f"{base_name}.svg")
    try:
        plot.save(svg_path, width=width, height=height, verbose=False)
        print(f"   Saved: {svg_path}")
    except Exception as exc:
        logger.debug("SVG export failed for %s: %s", base_name, exc)
        print("   (SVG export failed – PNG available)")
