"""
Parse differential region results from ChIP-Atlas Diff Analysis.

The primary output file (``wabi_result.bed``) is a tab-separated BED with
8 columns and no header::

    chrom   chromStart  chromEnd  counts_A  counts_B  logFC  pvalue  qvalue

``counts_A`` and ``counts_B`` are comma-separated normalised read counts
per experiment (e.g. ``"47.80,31.04,36.33"``).

The IGV BED file (``wabi_result.igv.bed``) is BED9+GFF3 for visualisation,
with itemRgb colours indicating direction:

    - Orange ``(222,131,68)``  → A-enriched (positive logFC)
    - Blue   ``(106,153,208)`` → B-enriched (negative logFC)
"""

from __future__ import annotations

import logging
from typing import List

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# Canonical column order for the 8-column result BED
RESULT_COLUMNS: List[str] = [
    "chrom",
    "chromStart",
    "chromEnd",
    "counts_a",
    "counts_b",
    "logFC",
    "pvalue",
    "qvalue",
]


# ─────────────────────────── Public API ────────────────────────────────────── #


def parse_bed_results(bed_path: str) -> pd.DataFrame:
    """
    Parse differential regions from the primary ChIP-Atlas BED output.

    Parameters
    ----------
    bed_path : str
        Path to ``wabi_result.bed`` from a completed job.

    Returns
    -------
    pd.DataFrame
        Columns:

        - ``chrom``, ``chromStart``, ``chromEnd`` – genomic coordinates
        - ``counts_a``, ``counts_b`` – raw comma-separated count strings
        - ``logFC`` – log₂ fold change (positive → A enriched)
        - ``pvalue`` – edgeR / metilene p-value
        - ``qvalue`` – Benjamini-Hochberg FDR
        - ``region_size`` – ``chromEnd − chromStart`` (bp)
        - ``direction`` – ``'A_enriched'`` or ``'B_enriched'``
        - ``significant`` – ``qvalue < 0.05``
        - ``mean_count_a``, ``mean_count_b`` – mean normalised counts
        - ``score`` – ``−log10(qvalue) × 100``, capped at 1 000 (BED-compatible)

    Raises
    ------
    RuntimeError
        If the file cannot be parsed.
    """
    try:
        df = pd.read_csv(
            bed_path,
            sep="\t",
            header=None,
            comment="#",
            dtype=str,
        )
    except pd.errors.EmptyDataError as exc:
        raise RuntimeError(f"BED file is empty: {bed_path}") from exc
    except pd.errors.ParserError as exc:
        raise RuntimeError(f"Failed to parse BED file '{bed_path}': {exc}") from exc

    if df.empty:
        raise RuntimeError(f"BED file is empty: {bed_path}")

    # Drop browser/track header lines
    mask = ~df.iloc[:, 0].str.startswith(("browser", "track"), na=False)
    df = df[mask].reset_index(drop=True)

    n_cols = len(df.columns)
    if n_cols == 8:
        df.columns = RESULT_COLUMNS
    elif n_cols == 9:
        # IGV BED9 format – parse differently
        df.columns = [
            "chrom", "chromStart", "chromEnd", "name", "score",
            "strand", "thickStart", "thickEnd", "itemRgb",
        ]
        return _parse_igv_bed(df)
    else:
        cols = RESULT_COLUMNS[: min(n_cols, len(RESULT_COLUMNS))]
        if n_cols > len(RESULT_COLUMNS):
            cols += [f"extra_{i}" for i in range(n_cols - len(RESULT_COLUMNS))]
        df.columns = cols

    # Numeric conversion
    for col in ("chromStart", "chromEnd", "logFC", "pvalue", "qvalue"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # Derived columns
    df["region_size"] = df["chromEnd"] - df["chromStart"]
    df["direction"]   = np.where(df["logFC"] > 0, "A_enriched", "B_enriched")
    df["significant"] = df["qvalue"] < 0.05
    df["mean_count_a"] = df["counts_a"].apply(_mean_counts)
    df["mean_count_b"] = df["counts_b"].apply(_mean_counts)
    df["score"] = np.minimum(
        -np.log10(df["qvalue"].clip(lower=1e-300)) * 100, 1000
    ).fillna(0).astype(int)

    # Sort by significance then effect size
    df["_abs_logFC"] = df["logFC"].abs()
    df = df.sort_values(
        ["qvalue", "_abs_logFC", "chrom", "chromStart"],
        ascending=[True, False, True, True],
    ).reset_index(drop=True)
    df = df.drop(columns=["_abs_logFC"])

    logger.info(
        "Parsed %d regions (%d significant at FDR < 0.05) from %s",
        len(df),
        int(df["significant"].sum()),
        bed_path,
    )
    return df


def summarize_regions(df: pd.DataFrame) -> None:
    """
    Print a concise summary of parsed differential regions to stdout.

    Parameters
    ----------
    df : pd.DataFrame
        Output of :func:`parse_bed_results`.
    """
    n_total = len(df)
    n_a     = int((df["direction"] == "A_enriched").sum()) if n_total else 0
    n_b     = int((df["direction"] == "B_enriched").sum()) if n_total else 0
    n_sig   = int(df["significant"].sum())                 if n_total else 0

    print(f"\n   Differential regions summary:")
    print(f"   Total:                  {n_total:,}")
    print(f"   A-enriched:             {n_a:,}")
    print(f"   B-enriched:             {n_b:,}")
    print(f"   Significant (FDR<0.05): {n_sig:,}")

    if n_sig and "direction" in df.columns:
        sig = df[df["significant"]]
        print(f"     Sig. A-enriched:     {int((sig['direction']=='A_enriched').sum()):,}")
        print(f"     Sig. B-enriched:     {int((sig['direction']=='B_enriched').sum()):,}")

    if n_total and "region_size" in df.columns:
        print(f"   Median region size:     {df['region_size'].median():.0f} bp")
    if n_total and "logFC" in df.columns:
        print(f"   logFC range:            {df['logFC'].min():.2f} – {df['logFC'].max():.2f}")
    if n_total and "qvalue" in df.columns:
        print(f"   Min Q-value:            {df['qvalue'].min():.2e}")
    if n_total:
        top = df["chrom"].value_counts().head(5)
        print(f"   Top chromosomes:        {', '.join(f'{c}({n})' for c,n in top.items())}")


# ─────────────────────────── Private helpers ───────────────────────────────── #


def _parse_igv_bed(df: pd.DataFrame) -> pd.DataFrame:
    """Parse IGV BED9+GFF3 as a fallback when the main BED is unavailable."""
    df["chromStart"] = pd.to_numeric(df["chromStart"], errors="coerce")
    df["chromEnd"]   = pd.to_numeric(df["chromEnd"],   errors="coerce")
    df["region_size"] = df["chromEnd"] - df["chromStart"]

    if "itemRgb" in df.columns:
        df["direction"] = df["itemRgb"].apply(_direction_from_rgb)
    else:
        df["direction"] = "unknown"

    df["logFC"]  = df["name"].apply(lambda x: _extract_gff3_field(x, "LogFC"))
    df["pvalue"] = df["name"].apply(lambda x: _extract_gff3_field(x, "P-value"))
    df["qvalue"] = df["name"].apply(lambda x: _extract_gff3_field(x, "Q-value"))
    df["significant"] = df["qvalue"] < 0.05
    df["score"] = pd.to_numeric(df.get("score", 0), errors="coerce").fillna(0)
    return df


def _extract_gff3_field(name_str: str, key: str) -> float:
    """Extract a numeric value from a URL-encoded GFF3 name field."""
    from urllib.parse import unquote
    try:
        decoded = unquote(str(name_str))
        for part in decoded.split(";"):
            if part.startswith(f"{key}="):
                return float(part.split("=", 1)[1])
    except (ValueError, IndexError):
        pass
    return np.nan


def _direction_from_rgb(item_rgb: str) -> str:
    """Infer enrichment direction from itemRgb colour string."""
    try:
        r, _g, b = (int(x) for x in str(item_rgb).split(","))
        if r > b:
            return "A_enriched"
        if b > r:
            return "B_enriched"
        return "unchanged"
    except (ValueError, IndexError):
        return "unknown"


def _mean_counts(counts_str: str) -> float:
    """Return mean of a comma-separated count string."""
    try:
        vals = [float(x) for x in str(counts_str).split(",")]
        return sum(vals) / len(vals) if vals else 0.0
    except (ValueError, TypeError):
        return 0.0
