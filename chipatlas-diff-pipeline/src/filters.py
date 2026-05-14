"""
Post-retrieval filtering of ChIP-Atlas differential regions.

Provides both automated QC filters (non-standard contigs, tiny regions)
and flexible user-facing filters (FDR, logFC, direction, chromosome, size).
"""

from __future__ import annotations

import logging
from typing import List, Optional

import pandas as pd

from src.validation import get_standard_chroms

logger = logging.getLogger(__name__)


# ─────────────────────────── QC filters ────────────────────────────────────── #


def filter_standard_chromosomes(
    df: pd.DataFrame,
    genome: Optional[str] = None,
) -> pd.DataFrame:
    """
    Remove regions on non-standard contigs (random scaffolds, unplaced seqs).

    Parameters
    ----------
    df : pd.DataFrame
        Parsed regions DataFrame.
    genome : str or None
        Assembly name. Falls back to hg38 chromosome set if None.

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame, index reset.
    """
    standard = get_standard_chroms(genome)
    before = len(df)
    filtered = df[df["chrom"].isin(standard)].copy()
    removed = before - len(filtered)
    if removed:
        logger.info("Removed %d regions on non-standard contigs.", removed)
        print(f"   Removed {removed:,} regions on non-standard contigs")
    return filtered.reset_index(drop=True)


def filter_min_region_size(
    df: pd.DataFrame,
    min_size: int = 10,
) -> pd.DataFrame:
    """
    Remove biologically implausible sub-10-bp regions.

    Parameters
    ----------
    df : pd.DataFrame
        Parsed regions DataFrame.
    min_size : int
        Minimum region length in bp (default 10).

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame, index reset.
    """
    if "region_size" not in df.columns:
        return df
    before = len(df)
    filtered = df[df["region_size"] >= min_size].copy()
    removed = before - len(filtered)
    if removed:
        logger.info("Removed %d regions < %d bp.", removed, min_size)
        print(f"   Removed {removed:,} regions < {min_size} bp")
    return filtered.reset_index(drop=True)


# ─────────────────────────── User filters ──────────────────────────────────── #


def filter_regions(
    regions_df: pd.DataFrame,
    max_qvalue: Optional[float] = None,
    min_logfc: Optional[float] = None,
    min_score: int = 0,
    min_region_size: int = 0,
    max_region_size: Optional[int] = None,
    direction_filter: Optional[str] = None,
    chromosomes: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Apply user-defined post-hoc filters to a differential regions DataFrame.

    Parameters
    ----------
    regions_df : pd.DataFrame
        Parsed DataFrame from :func:`~src.parser.parse_bed_results`.
    max_qvalue : float or None
        Maximum FDR threshold (e.g. ``0.05``). None → no filter.
    min_logfc : float or None
        Minimum absolute logFC (e.g. ``1.0``). None → no filter.
    min_score : int
        Minimum BED score (0–1000, default 0).
    min_region_size : int
        Minimum region size in bp (default 0).
    max_region_size : int or None
        Maximum region size in bp. None → no limit.
    direction_filter : str or None
        ``'A_enriched'`` or ``'B_enriched'``. None → keep both.
    chromosomes : list of str or None
        Whitelist of chromosomes (e.g. ``['chr1', 'chr2']``). None → all.

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame, index reset.
    """
    df = regions_df.copy()
    initial = len(df)

    if max_qvalue is not None and "qvalue" in df.columns:
        df = df[df["qvalue"] < max_qvalue]
        print(f"   After FDR < {max_qvalue}: {len(df):,} regions")

    if min_logfc is not None and "logFC" in df.columns:
        df = df[df["logFC"].abs() >= min_logfc]
        print(f"   After |logFC| ≥ {min_logfc}: {len(df):,} regions")

    if min_score > 0 and "score" in df.columns:
        df = df[df["score"] >= min_score]
        print(f"   After score ≥ {min_score}: {len(df):,} regions")

    if min_region_size > 0 and "region_size" in df.columns:
        df = df[df["region_size"] >= min_region_size]
        print(f"   After size ≥ {min_region_size} bp: {len(df):,} regions")

    if max_region_size is not None and "region_size" in df.columns:
        df = df[df["region_size"] <= max_region_size]
        print(f"   After size ≤ {max_region_size} bp: {len(df):,} regions")

    if direction_filter is not None and "direction" in df.columns:
        df = df[df["direction"] == direction_filter]
        print(f"   After direction = '{direction_filter}': {len(df):,} regions")

    if chromosomes and "chrom" in df.columns:
        df = df[df["chrom"].isin(chromosomes)]
        print(f"   After chrom filter ({chromosomes}): {len(df):,} regions")

    if len(df) < initial:
        logger.info("filter_regions: %d → %d regions.", initial, len(df))

    return df.reset_index(drop=True)
