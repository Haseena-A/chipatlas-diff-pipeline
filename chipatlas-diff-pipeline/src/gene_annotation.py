"""
Annotate differential regions with nearest/overlapping genes via UCSC REST API.

Queries the UCSC Genome Browser REST API to find genes overlapping or near
each differential region.  Up to *max_regions* (default 50) are annotated;
remaining rows receive an empty string.

Graceful fallback: if the UCSC API is unavailable, annotations are silently
skipped (empty strings returned) without raising exceptions.

UCSC REST API reference:
    https://genome.ucsc.edu/goldenPath/help/api.html
"""

from __future__ import annotations

import logging
from typing import Optional

import pandas as pd
import requests

logger = logging.getLogger(__name__)

UCSC_API_URL = "https://api.genome.ucsc.edu/getData/track"
FLANK_BP     = 5_000     # Search window around region when no direct overlap


def annotate_nearest_genes(
    df: pd.DataFrame,
    genome: str = "hg38",
    max_regions: int = 50,
) -> pd.Series:
    """
    Annotate genomic regions with the nearest gene symbol.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with ``chrom``, ``chromStart``, ``chromEnd`` columns.
        Should be pre-sorted by significance (most significant first).
    genome : str
        UCSC assembly name (``'hg38'``, ``'mm10'``, etc.).
    max_regions : int
        Maximum number of regions to annotate (default 50).

    Returns
    -------
    pd.Series
        Gene annotations indexed to match *df*.  Each entry is one of:

        * ``"GENE"``              – overlapping gene
        * ``"GENE (+2.1kb)"``    – nearest gene within *FLANK_BP*
        * ``""``                  – no nearby gene found / API unavailable
    """
    annotations = pd.Series("", index=df.index)
    subset = df.head(max_regions)
    if subset.empty:
        return annotations

    print(f"   Annotating top {len(subset)} regions with nearest genes (UCSC API) …")
    n_annotated = n_failed = 0

    for idx, row in subset.iterrows():
        try:
            gene_str = _query_region(
                chrom=str(row["chrom"]),
                start=int(row["chromStart"]),
                end=int(row["chromEnd"]),
                genome=genome,
            )
            if gene_str:
                annotations.at[idx] = gene_str
                n_annotated += 1
        except Exception as exc:
            n_failed += 1
            logger.debug("UCSC query failed for %s: %s", idx, exc)
            if n_failed >= 3:
                print(
                    f"   Warning: UCSC API unavailable ({n_failed} failures). "
                    "Skipping remaining gene annotations."
                )
                break

    print(f"   Annotated {n_annotated}/{len(subset)} regions with gene symbols")
    return annotations


# ─────────────────────────── Internal helpers ──────────────────────────────── #


def _query_region(chrom: str, start: int, end: int, genome: str) -> str:
    """
    Return a gene annotation string for one genomic region.

    Tries direct overlap first; falls back to nearest gene within FLANK_BP.
    """
    # Direct overlap
    genes = _fetch_gene_symbols(chrom, start, end, genome)
    if genes:
        return ", ".join(sorted(genes))

    # Nearest within flank
    exp_start = max(0, start - FLANK_BP)
    exp_end   = end + FLANK_BP
    records   = _fetch_gene_records(chrom, exp_start, exp_end, genome)
    if not records:
        return ""

    region_mid = (start + end) / 2
    best_gene:  Optional[str]   = None
    best_dist:  float           = float("inf")

    for rec in records:
        tx_start = rec.get("txStart", 0)
        tx_end   = rec.get("txEnd",   0)
        gene_mid = (tx_start + tx_end) / 2
        dist     = abs(region_mid - gene_mid)
        name2    = rec.get("name2", "")
        if name2 and dist < best_dist:
            best_dist = dist
            best_gene = name2

    if best_gene and best_dist > 0:
        return f"{best_gene} ({best_dist / 1000:+.1f}kb)"
    return best_gene or ""


def _fetch_gene_symbols(chrom: str, start: int, end: int, genome: str) -> set:
    """Return unique gene symbols overlapping a region."""
    records = _fetch_gene_records(chrom, start, end, genome)
    return {rec.get("name2", "") for rec in records if rec.get("name2")}


def _fetch_gene_records(chrom: str, start: int, end: int, genome: str) -> list:
    """Query the UCSC ncbiRefSeq track for a genomic interval."""
    url = (
        f"{UCSC_API_URL}?genome={genome};track=ncbiRefSeq;"
        f"chrom={chrom};start={int(start)};end={int(end)}"
    )
    try:
        resp = requests.get(url, timeout=10)
        resp.raise_for_status()
        return resp.json().get("ncbiRefSeq", [])
    except (requests.RequestException, ValueError, KeyError):
        return []
