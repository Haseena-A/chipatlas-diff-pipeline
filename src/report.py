"""
Markdown and summary report generation for ChIP-Atlas differential analysis.

Produces:
    - ``summary_report.md``  – human-readable Markdown report
    - ``diff_regions_all.csv``        – all differential regions (QC-filtered)
    - ``diff_regions_significant.csv`` – FDR < 0.05 regions
    - ``diff_regions_top50.csv``       – top 50 by significance (with gene annotations)
    - ``diff_regions_unfiltered.csv``  – pre-QC regions (if different)
    - ``analysis_object.pkl``          – complete results for downstream Python use
"""

from __future__ import annotations

import logging
import os
import pickle
from datetime import datetime
from typing import Dict, List, Optional

import pandas as pd

logger = logging.getLogger(__name__)


# ─────────────────────────── Entry point ───────────────────────────────────── #


def export_all(
    results: Dict,
    output_dir: str = "diff_analysis_results",
    qvalue_threshold: float = 0.05,
    annotate_genes: bool = True,
) -> None:
    """
    Export all results: CSVs, pickle, and Markdown report.

    Parameters
    ----------
    results : dict
        Output of :func:`~src.diff_analysis.run_diff_workflow`.
    output_dir : str
        Destination directory.
    qvalue_threshold : float
        FDR threshold for the *significant* subset (default 0.05).
    annotate_genes : bool
        Attempt UCSC gene annotation for top-50 regions (default True).
    """
    os.makedirs(output_dir, exist_ok=True)
    print("\n" + "=" * 70)
    print("EXPORTING RESULTS")
    print("=" * 70 + "\n")

    df           = results["diff_regions"]
    params       = results["parameters"]
    unfiltered   = results.get("diff_regions_unfiltered", df)
    qc_warnings  = results.get("qc_warnings", [])

    desc_a        = params.get("description_a", "Group A")
    desc_b        = params.get("description_b", "Group B")
    analysis_type = params.get("analysis_type", "diffbind")
    genome        = params.get("genome", "hg38")

    n_a_enriched = int((df["direction"] == "A_enriched").sum()) if len(df) else 0
    n_b_enriched = int((df["direction"] == "B_enriched").sum()) if len(df) else 0

    # ── 1. Pickle ────────────────────────────────────────────────────────── #
    print("1. Saving analysis object for downstream use …")
    obj = {
        "diff_regions":           df,
        "diff_regions_unfiltered": unfiltered,
        "qc_warnings":            qc_warnings,
        "experiments_a":          results["experiments_a"],
        "experiments_b":          results["experiments_b"],
        "raw_files":              results.get("raw_files", {}),
        "log_content":            results.get("log_content", ""),
        "parameters":             params,
        "n_regions_total":        len(df),
        "n_regions_unfiltered":   len(unfiltered),
        "n_regions_a_enriched":   n_a_enriched,
        "n_regions_b_enriched":   n_b_enriched,
        "timestamp":              datetime.now().isoformat(),
    }
    pkl_path = os.path.join(output_dir, "analysis_object.pkl")
    with open(pkl_path, "wb") as fh:
        pickle.dump(obj, fh)
    print(f"   Saved: {pkl_path}")
    print("   (Load with: import pickle; obj = pickle.load(open('analysis_object.pkl','rb')))")

    # ── 2. All regions ───────────────────────────────────────────────────── #
    print("\n2. Exporting differential regions …")
    all_path = os.path.join(output_dir, "diff_regions_all.csv")
    df.to_csv(all_path, index=False)
    print(f"   Saved: {all_path} ({len(df):,} regions)")

    # ── 3. Significant regions ───────────────────────────────────────────── #
    if "qvalue" in df.columns:
        sig = df[df["qvalue"] < qvalue_threshold]
        if len(sig):
            sig_path = os.path.join(output_dir, "diff_regions_significant.csv")
            sig.to_csv(sig_path, index=False)
            print(f"   Saved: {sig_path} ({len(sig):,} FDR < {qvalue_threshold} regions)")
        else:
            print(f"   No regions with FDR < {qvalue_threshold}")
    elif "score" in df.columns:
        sig = df[df["score"] >= 200]
        if len(sig):
            sig_path = os.path.join(output_dir, "diff_regions_significant.csv")
            sig.to_csv(sig_path, index=False)
            print(f"   Saved: {sig_path} ({len(sig):,} high-score regions)")

    # ── 4. Top 50 (with optional gene annotation) ────────────────────────── #
    top50 = _get_top50(df)
    gene_annotations: Optional[pd.Series] = None

    if annotate_genes and len(top50):
        try:
            from src.gene_annotation import annotate_nearest_genes
            gene_annotations = annotate_nearest_genes(top50, genome=genome, max_regions=50)
        except Exception as exc:
            logger.warning("Gene annotation skipped: %s", exc)
            print(f"   Warning: Gene annotation failed ({exc}). Skipping.")

    if gene_annotations is not None and len(gene_annotations):
        top50 = top50.copy()
        top50["nearest_gene"] = gene_annotations.reindex(top50.index).fillna("")

    top50_path = os.path.join(output_dir, "diff_regions_top50.csv")
    top50.to_csv(top50_path, index=False)
    print(f"   Saved: {top50_path} (top {len(top50)} by significance)")

    # ── 4b. Unfiltered ───────────────────────────────────────────────────── #
    if len(unfiltered) != len(df):
        unfilt_path = os.path.join(output_dir, "diff_regions_unfiltered.csv")
        unfiltered.to_csv(unfilt_path, index=False)
        print(f"   Saved: {unfilt_path} ({len(unfiltered):,} pre-QC regions)")

    # ── 5. Markdown report ───────────────────────────────────────────────── #
    print("\n3. Generating summary report …")
    report_path = os.path.join(output_dir, "summary_report.md")
    _write_markdown_report(
        path=report_path,
        results=results,
        df=df,
        unfiltered=unfiltered,
        top50=top50,
        gene_annotations=gene_annotations,
        qc_warnings=qc_warnings,
        params=params,
        desc_a=desc_a,
        desc_b=desc_b,
        analysis_type=analysis_type,
        n_a_enriched=n_a_enriched,
        n_b_enriched=n_b_enriched,
        qvalue_threshold=qvalue_threshold,
    )
    print(f"   Saved: {report_path}")

    # ── Final summary ────────────────────────────────────────────────────── #
    print("\n" + "=" * 70)
    print("=== Export Complete ===")
    print("=" * 70)
    print(f"\nAll results saved to: {output_dir}/")
    print(f"  analysis_object.pkl   – complete results for downstream use")
    print(f"  diff_regions_all.csv  – {len(df):,} QC-filtered regions")
    if "qvalue" in df.columns:
        n_sig = int((df["qvalue"] < qvalue_threshold).sum())
        if n_sig:
            print(f"  diff_regions_significant.csv – {n_sig:,} FDR < {qvalue_threshold}")
    print(f"  diff_regions_top50.csv – top {len(top50)} by significance")
    print(f"  summary_report.md")
    print()


# ─────────────────────────── Helpers ───────────────────────────────────────── #


def _get_top50(df: pd.DataFrame) -> pd.DataFrame:
    """Return top-50 rows sorted by significance then effect size."""
    if len(df) == 0:
        return df
    if "qvalue" in df.columns and "logFC" in df.columns:
        tmp = df.copy()
        tmp["_abs"] = tmp["logFC"].abs()
        top = tmp.sort_values(
            ["qvalue", "_abs", "chrom", "chromStart"],
            ascending=[True, False, True, True],
        ).head(50).drop(columns=["_abs"])
    elif "score" in df.columns:
        top = df.sort_values("score", ascending=False).head(50)
    else:
        top = df.head(50)
    return top


def _write_markdown_report(
    path: str,
    results: Dict,
    df: pd.DataFrame,
    unfiltered: pd.DataFrame,
    top50: pd.DataFrame,
    gene_annotations: Optional[pd.Series],
    qc_warnings: List[Dict],
    params: Dict,
    desc_a: str,
    desc_b: str,
    analysis_type: str,
    n_a_enriched: int,
    n_b_enriched: int,
    qvalue_threshold: float,
) -> None:
    """Write the Markdown summary report to *path*."""
    type_label = (
        "Differential Peak Regions"
        if analysis_type == "diffbind"
        else "Differentially Methylated Regions"
    )
    n_a_exp = len(results["experiments_a"])
    n_b_exp = len(results["experiments_b"])

    with open(path, "w") as fh:

        fh.write(f"# ChIP-Atlas Differential Analysis Summary\n\n")
        fh.write(f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        # Parameters
        fh.write("## Analysis Parameters\n\n")
        fh.write(f"| Parameter | Value |\n|-----------|-------|\n")
        fh.write(f"| Analysis type | {type_label} |\n")
        fh.write(f"| Genome | {params.get('genome', 'N/A')} |\n")
        fh.write(
            f"| Group A ({desc_a}) | "
            f"{', '.join(results['experiments_a'])} (n={n_a_exp}) |\n"
        )
        fh.write(
            f"| Group B ({desc_b}) | "
            f"{', '.join(results['experiments_b'])} (n={n_b_exp}) |\n\n"
        )

        # Experimental design caveats
        design_caveats = params.get("design_caveats", [])
        if design_caveats:
            fh.write("## ⚠️ Experimental Design Caveats\n\n")
            for c in design_caveats:
                fh.write(f"- {c}\n")
            fh.write("\n")

        # QC warnings
        if qc_warnings:
            fh.write("## QC Warnings\n\n")
            for w in qc_warnings:
                icon = "**MAJOR**" if w["severity"] == "major" else "Minor"
                fh.write(f"### {icon}: {w['issue']}\n\n")
                fh.write(f"{w['details']}\n\n")
                fh.write(f"**Recommendation:** {w['recommendation']}\n\n")

        # QC filtering note
        if len(unfiltered) != len(df):
            fh.write("## QC Filtering Applied\n\n")
            fh.write(f"- **Before QC filtering:** {len(unfiltered):,} regions\n")
            fh.write(f"- **After QC filtering:** {len(df):,} regions\n")
            fh.write(
                f"- **Removed:** {len(unfiltered) - len(df):,} regions "
                "(non-standard contigs and regions < 10 bp)\n"
            )
            fh.write("- **Unfiltered data:** `diff_regions_unfiltered.csv`\n\n")

        # Statistical caveats
        fh.write("## Statistical Caveats\n\n")
        fh.write(f"- **Sample size:** n={n_a_exp} (Group A) vs n={n_b_exp} (Group B)\n")
        if min(n_a_exp, n_b_exp) < 5:
            fh.write(
                "- **Power limitation:** With n < 5 per group, edgeR FDR estimates "
                "may be unstable for moderate effect sizes\n"
            )
        fh.write(
            "- **Validation:** Key findings should be validated with orthogonal "
            "methods or additional replicates\n"
        )
        fh.write("- **Public data:** Results depend on the quality of public ChIP-Atlas data\n\n")

        # Summary statistics
        fh.write("## Summary Statistics\n\n")
        fh.write(f"- **Total differential regions:** {len(df):,}\n")
        fh.write(f"- **Enriched in {desc_a}:** {n_a_enriched:,}\n")
        fh.write(f"- **Enriched in {desc_b}:** {n_b_enriched:,}\n")

        if "significant" in df.columns and len(df):
            n_sig = int(df["significant"].sum())
            fh.write(f"- **Significant (FDR < 0.05):** {n_sig:,}\n")
            if "direction" in df.columns and n_sig:
                sig_df = df[df["significant"]]
                fh.write(
                    f"  - {desc_a}-enriched: "
                    f"{int((sig_df['direction']=='A_enriched').sum()):,}\n"
                )
                fh.write(
                    f"  - {desc_b}-enriched: "
                    f"{int((sig_df['direction']=='B_enriched').sum()):,}\n"
                )
        if "logFC" in df.columns and len(df):
            fh.write(
                f"- **logFC range:** {df['logFC'].min():.2f} to {df['logFC'].max():.2f}\n"
            )
        if "qvalue" in df.columns and len(df):
            fh.write(f"- **Min Q-value:** {df['qvalue'].min():.2e}\n")
        if "region_size" in df.columns and len(df):
            fh.write(f"- **Median region size:** {df['region_size'].median():.0f} bp\n")
        fh.write("\n")

        # Top 10 table
        if len(df):
            has_genes = gene_annotations is not None and len(gene_annotations) > 0
            gene_lookup = (
                {i: g for i, g in gene_annotations.items() if g}
                if has_genes else {}
            )
            top10 = _get_top50(df).head(10)

            fh.write("## Top 10 Differential Regions\n\n")
            if has_genes:
                fh.write(
                    "| Rank | Location | Size (bp) | logFC | Q-value | "
                    "Direction | Nearest Gene |\n"
                )
                fh.write(
                    "|------|----------|-----------|-------|---------|"
                    "----------|--------------|\n"
                )
            else:
                fh.write(
                    "| Rank | Location | Size (bp) | logFC | Q-value | Direction |\n"
                )
                fh.write(
                    "|------|----------|-----------|-------|---------|----------|\n"
                )

            for rank, (idx, row) in enumerate(top10.iterrows(), 1):
                loc  = (
                    f"{row['chrom']}:{int(row['chromStart']):,}-{int(row['chromEnd']):,}"
                )
                size = int(row.get("region_size", row["chromEnd"] - row["chromStart"]))
                lfc  = f"{row['logFC']:.2f}"  if "logFC"  in row else "N/A"
                qval = f"{row['qvalue']:.2e}" if "qvalue" in row else "N/A"
                dirn = row.get("direction", "unknown")
                if has_genes:
                    gene = gene_lookup.get(idx, "")
                    fh.write(
                        f"| {rank} | {loc} | {size:,} | {lfc} | {qval} "
                        f"| {dirn} | {gene} |\n"
                    )
                else:
                    fh.write(
                        f"| {rank} | {loc} | {size:,} | {lfc} | {qval} | {dirn} |\n"
                    )
            fh.write("\n")

        # Chromosome distribution
        if len(df):
            fh.write("## Chromosome Distribution\n\n")
            for chrom, count in df["chrom"].value_counts().head(10).items():
                fh.write(f"- **{chrom}:** {count:,} regions\n")
            fh.write("\n")

        # Files generated
        fh.write("## Files Generated\n\n")
        fh.write("**Analysis object:**\n")
        fh.write("- `analysis_object.pkl` – complete results for downstream Python use\n\n")
        fh.write("**CSV tables:**\n")
        fh.write("- `diff_regions_all.csv` – all QC-filtered differential regions\n")
        fh.write(f"- `diff_regions_significant.csv` – FDR < {qvalue_threshold} regions\n")
        fh.write("- `diff_regions_top50.csv` – top 50 by significance\n")
        if len(unfiltered) != len(df):
            fh.write("- `diff_regions_unfiltered.csv` – pre-QC filter regions\n")
        fh.write("\n**Visualisations (PNG + SVG, 300 DPI):**\n")
        fh.write("- `volcano_plot` – logFC vs −log₁₀(Q-value)\n")
        fh.write("- `chromosome_distribution` – significant regions by chromosome\n")
        fh.write("- `region_size_distribution` – region size histogram\n")
        fh.write("- `ma_plot` – log₂ mean count vs logFC\n\n")
        fh.write("**Raw results:**\n")
        fh.write("- `raw_results/` – original BED, IGV BED, log, and XML files\n\n")

        # Next steps
        fh.write("## Suggested Next Steps\n\n")
        fh.write(
            "1. **Visualise in IGV** – Open the `.igv.bed` file to inspect "
            "differential regions in their genomic context\n"
        )
        fh.write(
            "2. **Motif analysis** – Run HOMER or MEME-ChIP on significant regions "
            "to identify enriched TF motifs\n"
        )
        fh.write(
            "3. **Gene ontology** – Use GREAT or nearby gene annotations to identify "
            "enriched biological processes\n"
        )
        fh.write(
            "4. **Integrate with expression data** – Cross-reference peaks with "
            "RNA-seq DE results to link regulatory regions to gene expression\n"
        )
