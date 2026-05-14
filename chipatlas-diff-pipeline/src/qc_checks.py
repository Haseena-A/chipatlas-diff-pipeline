"""
Automated QC checks for ChIP-Atlas differential analysis results.

Flags interpretive pitfalls that commonly arise with public ChIP-seq data:

* Sex-chromosome confounds (chrY/chrX from mismatched sample sex)
* Non-standard contigs (mapping artifacts on unplaced scaffolds)
* Implausibly small regions (< 10 bp — likely edge artifacts)
* Mitochondrial DNA enrichment (non-specific signal for ChIP/ATAC)
* Dense genomic clusters (possible CNV or structural-variant artefacts)
* MA-plot asymmetry (low-count noise driving one direction)
* Low sample size (limited statistical power)

Each check returns a list of warning dicts with keys:
    ``severity`` ('major' or 'minor'), ``issue``, ``details``,
    ``recommendation``, and optionally ``flagged_regions``.
"""

from __future__ import annotations

import logging
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from src.validation import get_standard_chroms

logger = logging.getLogger(__name__)

SEX_CHROMS = {"chrX", "chrY"}


# ─────────────────────────── Entry point ───────────────────────────────────── #


def run_qc_checks(
    df: pd.DataFrame,
    n_samples_a: int,
    n_samples_b: int,
    analysis_type: str = "diffbind",
    genome: Optional[str] = None,
) -> List[Dict]:
    """
    Run the full QC suite on parsed (unfiltered) differential regions.

    Parameters
    ----------
    df : pd.DataFrame
        Parsed DataFrame from :func:`~src.parser.parse_bed_results`.
    n_samples_a : int
        Number of experiments in group A.
    n_samples_b : int
        Number of experiments in group B.
    analysis_type : str
        ``'diffbind'`` or ``'dmr'``.
    genome : str or None
        Genome assembly for chromosome validation.

    Returns
    -------
    list of dict
        QC warnings (empty list = no issues detected).
    """
    warnings: List[Dict] = []
    if len(df) == 0:
        return warnings

    print("\n   Running QC checks …")

    warnings.extend(_check_sex_chromosome_confound(df))
    warnings.extend(_check_sex_chromosome_significant(df))
    warnings.extend(_check_nonstandard_contigs(df, genome=genome))
    warnings.extend(_check_small_regions(df))
    warnings.extend(_check_chrM_enrichment(df, analysis_type))
    warnings.extend(_check_genomic_clusters(df))
    warnings.extend(_check_ma_asymmetry(df))
    warnings.extend(_check_region_size_distribution(df, analysis_type))
    warnings.extend(_check_sample_size(n_samples_a, n_samples_b))

    # Deduplicate: prefer the "significant regions" variant when both exist
    sig_warned_chroms: set = set()
    for w in warnings:
        if "significant regions" in w.get("issue", ""):
            sig_warned_chroms.add(w["issue"].split(" ")[0])

    if sig_warned_chroms:
        warnings = [
            w for w in warnings
            if not (
                "sex chromosome confound" in w.get("issue", "")
                and "significant regions" not in w.get("issue", "")
                and w["issue"].split(" ")[0] in sig_warned_chroms
            )
        ]

    if warnings:
        print(f"\n   ⚠️  {len(warnings)} QC warning(s) detected:")
        for w in warnings:
            icon = "🔴" if w["severity"] == "major" else "🟡"
            print(f"   {icon} [{w['severity'].upper()}] {w['issue']}")
            print(f"      {w['details'][:120]}…")
            print(f"      → {w['recommendation'][:120]}")
    else:
        print("   ✓ No QC issues detected")

    return warnings


# ─────────────────────────── Individual checks ─────────────────────────────── #


def _check_sex_chromosome_confound(df: pd.DataFrame) -> List[Dict]:
    """Flag probable sex-chromosome artefacts across all regions."""
    warnings: List[Dict] = []
    for chrom in ("chrY", "chrX"):
        sub = df[df["chrom"] == chrom]
        if len(sub) < 3:
            continue
        n_a = int((sub["direction"] == "A_enriched").sum())
        n_b = int((sub["direction"] == "B_enriched").sum())
        n   = len(sub)
        dominant_frac = max(n_a, n_b) / n
        if dominant_frac < 0.8:
            continue
        if "mean_count_a" not in sub.columns or "mean_count_b" not in sub.columns:
            continue
        enriched_grp = "A" if n_a > n_b else "B"
        dep_col = "mean_count_b" if enriched_grp == "A" else "mean_count_a"
        enr_col = "mean_count_a" if enriched_grp == "A" else "mean_count_b"
        dep_mean = sub[dep_col].mean()
        enr_mean = sub[enr_col].mean()
        if enr_mean == 0 or dep_mean / enr_mean >= 0.1:
            continue
        n_sig = int(sub["significant"].sum()) if "significant" in sub.columns else 0
        warnings.append({
            "severity": "major",
            "issue": f"{chrom} sex chromosome confound",
            "details": (
                f"{n} {chrom} regions are {dominant_frac:.0%} enriched in group "
                f"{enriched_grp} with near-zero counts in the other group "
                f"(mean {dep_mean:.1f} vs {enr_mean:.1f}). "
                f"{n_sig} are statistically significant. "
                f"Pattern indicates sex mismatch, not genuine differential binding."
            ),
            "recommendation": (
                f"Exclude {chrom} regions from biological interpretation. "
                f"Report as probable sex-chromosome artefacts."
            ),
        })
    return warnings


def _check_sex_chromosome_significant(df: pd.DataFrame) -> List[Dict]:
    """Flag sex-chromosome artefacts among *significant* regions specifically."""
    warnings: List[Dict] = []
    if "significant" not in df.columns:
        return warnings
    for chrom in ("chrY", "chrX"):
        sig_sub = df[(df["chrom"] == chrom) & df["significant"]]
        if len(sig_sub) < 2:
            continue
        n_a = int((sig_sub["direction"] == "A_enriched").sum())
        n_b = int((sig_sub["direction"] == "B_enriched").sum())
        n_sig = len(sig_sub)
        dominant_frac = max(n_a, n_b) / n_sig
        if dominant_frac < 0.67:
            continue
        if "mean_count_a" not in sig_sub.columns:
            continue
        enriched_grp = "A" if n_a > n_b else "B"
        dep_col = "mean_count_b" if enriched_grp == "A" else "mean_count_a"
        enr_col = "mean_count_a" if enriched_grp == "A" else "mean_count_b"
        dom_dir = "A_enriched" if n_a > n_b else "B_enriched"
        dom_sig = sig_sub[sig_sub["direction"] == dom_dir]
        dep_mean = dom_sig[dep_col].mean()
        enr_mean = dom_sig[enr_col].mean()
        if enr_mean == 0 or dep_mean / enr_mean >= 0.1:
            continue
        warnings.append({
            "severity": "major",
            "issue": f"{chrom} sex chromosome confound (significant regions)",
            "details": (
                f"{max(n_a,n_b)}/{n_sig} significant {chrom} regions enriched in "
                f"group {enriched_grp} with near-zero counts in the other group "
                f"(mean {dep_mean:.1f} vs {enr_mean:.1f})."
            ),
            "recommendation": (
                f"Exclude significant {chrom} regions from biological interpretation."
            ),
        })
    return warnings


def _check_nonstandard_contigs(
    df: pd.DataFrame,
    genome: Optional[str] = None,
) -> List[Dict]:
    """Warn about regions on unplaced / random scaffolds."""
    standard = get_standard_chroms(genome)
    mask = ~df["chrom"].isin(standard)
    n_ns = int(mask.sum())
    if n_ns == 0:
        return []
    n_sig_ns = int((mask & df.get("significant", pd.Series(False, index=df.index))).sum())
    pct = n_ns / len(df) * 100
    n_contigs = df.loc[mask, "chrom"].nunique()
    return [{
        "severity": "minor",
        "issue": "Non-standard contigs present",
        "details": (
            f"{n_ns} regions ({pct:.1f}%) on {n_contigs} non-standard contig(s). "
            f"{n_sig_ns} are statistically significant."
        ),
        "recommendation": (
            "Default QC filtering removes non-standard contigs. "
            "Unfiltered data available in diff_regions_unfiltered.csv."
        ),
    }]


def _check_small_regions(df: pd.DataFrame) -> List[Dict]:
    """Flag sub-10-bp regions as biologically implausible."""
    if "region_size" not in df.columns:
        return []
    mask = df["region_size"] < 10
    n_small = int(mask.sum())
    if n_small == 0:
        return []
    n_sig_small = int(
        (mask & df.get("significant", pd.Series(False, index=df.index))).sum()
    )
    min_size = int(df.loc[mask, "region_size"].min())
    return [{
        "severity": "minor",
        "issue": "Implausibly small regions detected",
        "details": (
            f"{n_small} regions < 10 bp (smallest: {min_size} bp). "
            f"{n_sig_small} are statistically significant."
        ),
        "recommendation": (
            "Default QC filtering removes regions < 10 bp. "
            "Unfiltered data available in diff_regions_unfiltered.csv."
        ),
    }]


def _check_chrM_enrichment(
    df: pd.DataFrame,
    analysis_type: str = "diffbind",
) -> List[Dict]:
    """Flag chrM signal as probable contamination in ChIP/ATAC analyses."""
    if analysis_type == "dmr":
        return []
    sub = df[df["chrom"] == "chrM"]
    if len(sub) < 3:
        return []
    n_sig_chrm = int(sub["significant"].sum()) if "significant" in sub.columns else 0
    if n_sig_chrm == 0:
        return []
    n_sig_total = int(df["significant"].sum()) if "significant" in df.columns else 1
    pct = n_sig_chrm / n_sig_total * 100
    severity = "major" if (n_sig_chrm >= 5 or pct >= 10) else "minor"
    return [{
        "severity": severity,
        "issue": f"chrM enrichment ({n_sig_chrm} significant regions)",
        "details": (
            f"{n_sig_chrm} of {n_sig_total} significant regions ({pct:.0f}%) "
            f"are on chrM. Mitochondrial DNA lacks histones and is not a target "
            f"for most transcription factors — chrM ChIP signal typically reflects "
            f"non-specific antibody binding or library prep artefacts."
        ),
        "recommendation": (
            "Exclude chrM regions from biological interpretation. "
            "Consider filtering chrM from the analysis entirely."
        ),
    }]


def _check_genomic_clusters(
    df: pd.DataFrame,
    min_regions: int = 10,
    max_window: int = 2_000_000,
) -> List[Dict]:
    """Detect dense significant-region clusters (possible CNV artefacts)."""
    warnings: List[Dict] = []
    if "significant" not in df.columns:
        return warnings
    sig = df[df["significant"]].copy()
    if len(sig) < min_regions:
        return warnings

    for chrom in sig["chrom"].unique():
        csig = sig[sig["chrom"] == chrom].sort_values("chromStart")
        if len(csig) < min_regions:
            continue
        starts = csig["chromStart"].values
        ends   = csig["chromEnd"].values
        best_count, best_start, best_end = 0, 0, 0
        j = 0
        for i in range(len(starts)):
            while j < len(starts) and starts[j] - starts[i] <= max_window:
                j += 1
            count = j - i
            if count > best_count:
                best_count = count
                best_start = int(starts[i])
                best_end   = int(ends[j - 1])

        if best_count < min_regions:
            continue

        cluster = csig[
            (csig["chromStart"] >= best_start) &
            (csig["chromStart"] <= best_start + max_window)
        ]
        n_a = int((cluster["direction"] == "A_enriched").sum())
        n_b = int((cluster["direction"] == "B_enriched").sum())
        dom = "A" if n_a >= n_b else "B"
        dom_pct = max(n_a, n_b) / best_count * 100
        window_mb = (best_end - best_start) / 1e6
        flagged = list(zip(
            cluster["chrom"].values,
            cluster["chromStart"].values.astype(int),
            cluster["chromEnd"].values.astype(int),
        ))
        warnings.append({
            "severity": "minor",
            "issue": (
                f"{chrom} genomic cluster "
                f"({best_count} significant regions in {window_mb:.1f} Mb)"
            ),
            "details": (
                f"{best_count} significant regions cluster within "
                f"{chrom}:{best_start:,}-{best_end:,} ({window_mb:.1f} Mb). "
                f"{dom_pct:.0f}% enriched in group {dom} "
                f"({n_a} A-enriched, {n_b} B-enriched). "
                f"Dense clustering may indicate a copy number amplification "
                f"or structural variant rather than genuine differential binding."
            ),
            "recommendation": (
                f"Inspect {chrom}:{best_start:,}-{best_end:,} in a genome browser. "
                f"Check for known CNVs in these cell lines."
            ),
            "flagged_regions": flagged,
        })
    return warnings


def _check_ma_asymmetry(df: pd.DataFrame) -> List[Dict]:
    """Detect low-count noise driving one direction on the MA plot."""
    if "significant" not in df.columns:
        return []
    if "mean_count_a" not in df.columns or "mean_count_b" not in df.columns:
        return []
    sig = df[df["significant"]].copy()
    if len(sig) < 5:
        return []
    sig_a = sig[sig["direction"] == "A_enriched"]
    sig_b = sig[sig["direction"] == "B_enriched"]
    if len(sig_a) < 3 or len(sig_b) < 3:
        return []

    def _log2_mean(row: pd.Series) -> float:
        return float(np.log2(max((row["mean_count_a"] + row["mean_count_b"]) / 2, 0.1)))

    med_a = sig_a.apply(_log2_mean, axis=1).median()
    med_b = sig_b.apply(_log2_mean, axis=1).median()
    diff  = abs(med_a - med_b)
    if diff < 1.0:
        return []
    low_grp    = "A" if med_a < med_b else "B"
    low_med    = min(med_a, med_b)
    high_med   = max(med_a, med_b)
    low_n      = len(sig_a) if low_grp == "A" else len(sig_b)
    return [{
        "severity": "minor",
        "issue": f"MA plot asymmetry (group {low_grp} hits at low counts)",
        "details": (
            f"Group {low_grp} significant regions have a median log₂(mean count) "
            f"of {low_med:.1f} vs {high_med:.1f} for the other direction "
            f"(Δ={diff:.1f} log₂ units, {low_n} regions). "
            f"Asymmetric significant hits at low counts may reflect a "
            f"normalisation artefact rather than genuine differential binding."
        ),
        "recommendation": (
            f"Examine the MA plot. Focus biological interpretation on "
            f"higher-count significant regions in group {low_grp}."
        ),
    }]


def _check_region_size_distribution(
    df: pd.DataFrame,
    analysis_type: str = "diffbind",
) -> List[Dict]:
    """Warn when significant regions are unusually small (< 100 bp median)."""
    if analysis_type == "dmr":
        return []
    if "region_size" not in df.columns or "significant" not in df.columns:
        return []
    sig = df[df["significant"]]
    if len(sig) < 5:
        return []
    med = sig["region_size"].median()
    if med >= 100:
        return []
    return [{
        "severity": "minor",
        "issue": f"Unusually small significant regions (median {med:.0f} bp)",
        "details": (
            f"Median significant region size is {med:.0f} bp — smaller than "
            f"typical TF binding sites (100–500 bp) or histone marks (200 bp+). "
            f"May indicate peak fragmentation in the differential analysis."
        ),
        "recommendation": (
            "Examine representative regions in a genome browser. Sub-peaks may "
            "still be biologically relevant but should be interpreted cautiously."
        ),
    }]


def _check_sample_size(n_a: int, n_b: int) -> List[Dict]:
    """Warn about limited statistical power with small sample sizes."""
    min_n = min(n_a, n_b)
    if min_n <= 2:
        return [{
            "severity": "major",
            "issue": f"Very low sample size (n={n_a} vs n={n_b})",
            "details": (
                f"With n ≤ 3 per group, edgeR dispersion estimates are unreliable. "
                f"FDR values should be treated as approximate rankings."
            ),
            "recommendation": (
                "Focus on large effect sizes (|logFC| > 2). Validate with "
                "orthogonal methods or additional replicates."
            ),
        }]
    if min_n < 5:
        return [{
            "severity": "minor",
            "issue": f"Low sample size (n={n_a} vs n={n_b})",
            "details": (
                f"With n < 5 per group, edgeR FDR estimates can be unstable "
                f"for moderate effect sizes."
            ),
            "recommendation": (
                "Validate key findings with orthogonal methods or additional replicates."
            ),
        }]
    return []
