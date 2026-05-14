"""
End-to-end ChIP-Atlas differential analysis workflow orchestrator.

Chains together all pipeline stages:

    1. Input validation & pre-flight checks
    2. API job submission (with optional cache lookup)
    3. Status polling
    4. ZIP retrieval & extraction
    5. BED parsing
    6. QC checks (sex-chrom confounds, CNV clusters, low-count artefacts, …)
    7. Standard QC filtering (non-standard contigs, < 10 bp regions)
    8. User-requested post-filters
    9. Visualisation (volcano, chromosome, MA, size plots)
   10. Export (CSVs, pickle, Markdown report)
"""

from __future__ import annotations

import logging
import os
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)


def run_diff_workflow(
    experiments_a: List[str],
    experiments_b: List[str],
    genome: str = "hg38",
    analysis_type: str = "diffbind",
    title: str = "diff_analysis",
    description_a: str = "group_A",
    description_b: str = "group_B",
    min_score: int = 0,
    min_region_size: int = 0,
    direction_filter: Optional[str] = None,
    design_caveats: Optional[List[str]] = None,
    output_dir: str = "diff_analysis_results",
    use_cache: bool = True,
    generate_plots: bool = True,
    annotate_genes: bool = True,
    job_timeout: int = 900,
    poll_interval: int = 15,
) -> Dict:
    """
    Run the complete ChIP-Atlas differential analysis pipeline.

    Parameters
    ----------
    experiments_a : list of str
        SRX/ERX/DRX/GSM accessions for group A (minimum 2).
    experiments_b : list of str
        SRX/ERX/DRX/GSM accessions for group B (minimum 2).
    genome : str
        UCSC assembly (``'hg38'``, ``'mm10'``, …).
    analysis_type : str
        ``'diffbind'`` (edgeR DPR) or ``'dmr'`` (metilene DMR).
    title : str
        Short analysis title stored in job metadata.
    description_a : str
        Label for group A used in reports and plots.
    description_b : str
        Label for group B used in reports and plots.
    min_score : int
        Post-filter: minimum BED score (0–1000).
    min_region_size : int
        Post-filter: minimum region size in bp.
    direction_filter : str or None
        Post-filter: ``'A_enriched'`` or ``'B_enriched'``.
    design_caveats : list of str or None
        Free-text experimental design notes appended to the report.
    output_dir : str
        Root output directory.
    use_cache : bool
        Whether to check the local cache before submitting a new job.
    generate_plots : bool
        Whether to generate visualisations.
    annotate_genes : bool
        Whether to annotate the top-50 regions with UCSC gene symbols.
    job_timeout : int
        Maximum seconds to wait for the job to finish.
    poll_interval : int
        Seconds between status polls.

    Returns
    -------
    dict with keys
        - ``diff_regions`` – QC-filtered pd.DataFrame
        - ``diff_regions_unfiltered`` – pre-QC pd.DataFrame
        - ``qc_warnings`` – list of warning dicts
        - ``experiments_a``, ``experiments_b`` – input accession lists
        - ``raw_files`` – dict of extracted file paths
        - ``log_content`` – edgeR / metilene log string
        - ``parameters`` – echo of all run parameters
    """
    from src.api_client import submit_diff_job, poll_job_status, retrieve_zip_results
    from src.parser import parse_bed_results, summarize_regions
    from src.filters import filter_standard_chromosomes, filter_min_region_size, filter_regions
    from src.qc_checks import run_qc_checks
    from src.validation import run_preflight_checks
    from src.cache import AnalysisCache, make_cache_key
    from src.report import export_all
    from src.visualization import generate_all_plots

    # ── Validate ─────────────────────────────────────────────────────────── #
    run_preflight_checks(experiments_a, experiments_b, genome, analysis_type)
    os.makedirs(output_dir, exist_ok=True)

    analysis_labels = {
        "diffbind": "Differential Peak Regions (edgeR DPR)",
        "dmr":      "Differentially Methylated Regions (metilene)",
    }

    print("\n" + "=" * 70)
    print("CHIP-ATLAS DIFFERENTIAL ANALYSIS PIPELINE")
    print("=" * 70 + "\n")
    print(f"  Analysis type : {analysis_labels.get(analysis_type, analysis_type)}")
    print(f"  Genome        : {genome}")
    print(f"  Group A ({description_a}): {len(experiments_a)} experiment(s)")
    for eid in experiments_a:
        print(f"      - {eid}")
    print(f"  Group B ({description_b}): {len(experiments_b)} experiment(s)")
    for eid in experiments_b:
        print(f"      - {eid}")

    # ── Cache lookup ─────────────────────────────────────────────────────── #
    cache     = AnalysisCache()
    cache_key = make_cache_key(experiments_a, experiments_b, genome, analysis_type)

    if use_cache:
        cached = cache.get(cache_key)
        if cached is not None:
            print(f"\n✓ Cache hit ({cache_key}) – skipping API submission.")
            print("  (Delete ~/.chipatlas_cache/ to force re-run)\n")
            return cached

    # ── Step 1: Submit ───────────────────────────────────────────────────── #
    print("\nStep 1: Submitting job to ChIP-Atlas API …")
    request_id = submit_diff_job(
        experiments_a=experiments_a,
        experiments_b=experiments_b,
        genome=genome,
        analysis_type=analysis_type,
        title=title,
        description_a=description_a,
        description_b=description_b,
    )
    print(f"  Request ID: {request_id}")
    print("✓ Job submitted")

    # ── Step 2: Poll ─────────────────────────────────────────────────────── #
    print("\nStep 2: Waiting for analysis to complete …")
    print("  (Diff analysis typically takes 2–10 min depending on dataset size)")
    status = poll_job_status(
        request_id,
        timeout=job_timeout,
        poll_interval=poll_interval,
    )
    print(f"✓ Job finished (status: {status})")

    # ── Step 3: Retrieve & parse ─────────────────────────────────────────── #
    print("\nStep 3: Retrieving and parsing results …")
    extracted = retrieve_zip_results(request_id, output_dir)
    print(f"  Extracted {len(extracted)} file(s):")
    for k, p in extracted.items():
        print(f"    {k}: {os.path.basename(p)} ({os.path.getsize(p):,} bytes)")

    log_content = ""
    if "log" in extracted:
        with open(extracted["log"]) as fh:
            log_content = fh.read()

    results_df = parse_bed_results(extracted["bed"])
    n_raw = len(results_df)
    n_sig_raw = int(results_df["significant"].sum()) if "significant" in results_df.columns else 0
    print(f"✓ Parsed {n_raw:,} regions ({n_sig_raw:,} significant at FDR < 0.05)")

    # ── QC checks ────────────────────────────────────────────────────────── #
    qc_warnings = run_qc_checks(
        results_df,
        n_samples_a=len(experiments_a),
        n_samples_b=len(experiments_b),
        analysis_type=analysis_type,
        genome=genome,
    )

    # ── Standard QC filtering ─────────────────────────────────────────────── #
    unfiltered_df = results_df.copy()
    print("\n  Applying default QC filters …")
    results_df = filter_standard_chromosomes(results_df, genome=genome)
    results_df = filter_min_region_size(results_df, min_size=10)
    n_after_qc = len(results_df)
    if n_raw != n_after_qc:
        n_sig_after = int(results_df["significant"].sum()) if "significant" in results_df.columns else 0
        print(
            f"  QC filtering: {n_raw:,} → {n_after_qc:,} regions "
            f"({n_sig_raw - n_sig_after:,} significant removed)"
        )

    # ── Step 4: User post-filters ─────────────────────────────────────────── #
    if min_score or min_region_size or direction_filter:
        print("\nStep 4: Applying additional user filters …")
        results_df = filter_regions(
            results_df,
            min_score=min_score,
            min_region_size=min_region_size,
            direction_filter=direction_filter,
        )
    else:
        print("\nStep 4: No additional user filters requested")

    # ── Post-QC summary ───────────────────────────────────────────────────── #
    summarize_regions(results_df)

    n_final        = len(results_df)
    n_a_enriched   = int((results_df["direction"] == "A_enriched").sum()) if n_final else 0
    n_b_enriched   = int((results_df["direction"] == "B_enriched").sum()) if n_final else 0
    n_sig_final    = int(results_df["significant"].sum()) if n_final else 0

    print(f"\n  === Results Summary (post-QC) ===")
    print(f"  Total differential regions : {n_final:,}")
    print(f"  Significant (FDR < 0.05)   : {n_sig_final:,}")
    print(f"  Enriched in {description_a} : {n_a_enriched:,}")
    print(f"  Enriched in {description_b} : {n_b_enriched:,}")
    if n_final and "region_size" in results_df.columns:
        print(f"  Median region size         : {results_df['region_size'].median():.0f} bp")

    print("\n✓ Diff analysis completed successfully!")
    print("=" * 70 + "\n")

    # ── Assemble result dict ──────────────────────────────────────────────── #
    results = {
        "diff_regions":            results_df,
        "diff_regions_unfiltered": unfiltered_df,
        "qc_warnings":             qc_warnings,
        "experiments_a":           experiments_a,
        "experiments_b":           experiments_b,
        "raw_files":               extracted,
        "log_content":             log_content,
        "parameters": {
            "genome":           genome,
            "analysis_type":    analysis_type,
            "title":            title,
            "description_a":    description_a,
            "description_b":    description_b,
            "min_score":        min_score,
            "min_region_size":  min_region_size,
            "direction_filter": direction_filter,
            "design_caveats":   design_caveats or [],
        },
    }

    # ── Visualisations ────────────────────────────────────────────────────── #
    if generate_plots:
        generate_all_plots(results, output_dir=output_dir)

    # ── Export ────────────────────────────────────────────────────────────── #
    export_all(results, output_dir=output_dir, annotate_genes=annotate_genes)

    # ── Cache ─────────────────────────────────────────────────────────────── #
    if use_cache:
        cache.set(cache_key, results)

    return results
