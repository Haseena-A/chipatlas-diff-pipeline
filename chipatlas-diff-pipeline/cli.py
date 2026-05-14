#!/usr/bin/env python3
"""
chipatlas-diff  –  Command-line interface for the ChIP-Atlas differential
                   binding / methylation analysis pipeline.

Usage examples
--------------
# Minimal – two text files, one ID per line
chipatlas-diff run \\
    --group-a examples/group_a.txt \\
    --group-b examples/group_b.txt \\
    --genome hg38 \\
    --analysis-type diffbind \\
    --output-dir results/my_analysis

# YAML config
chipatlas-diff run --config configs/example_diffbind.yaml

# From example dataset (built-in)
chipatlas-diff example --dataset tp53_k562_dmso_vs_daunorubicin

# List cached analyses
chipatlas-diff cache list

# Clear cache
chipatlas-diff cache clear
"""

from __future__ import annotations

import argparse
import sys
import os


# ─────────────────────────── Sub-command handlers ─────────────────────────── #


def cmd_run(args: argparse.Namespace) -> int:
    """Handle the ``run`` sub-command."""
    import yaml  # optional dependency — only needed for --config

    from src.utils import setup_logging, load_ids_from_file, load_two_group_file
    from src.diff_analysis import run_diff_workflow

    setup_logging(level=args.log_level)

    # ── Load parameters ───────────────────────────────────────────────────── #
    params: dict = {}

    if args.config:
        with open(args.config) as fh:
            params = yaml.safe_load(fh)

    # CLI flags override config file values
    def _get(key, default=None):
        return getattr(args, key, None) or params.get(key, default)

    genome        = _get("genome", "hg38")
    analysis_type = _get("analysis_type", "diffbind")
    title         = _get("title", "chipatlas_diff")
    description_a = _get("description_a", "group_A")
    description_b = _get("description_b", "group_B")
    output_dir    = _get("output_dir", "diff_analysis_results")
    min_score     = int(_get("min_score", 0))
    min_region_size = int(_get("min_region_size", 0))
    direction_filter = _get("direction_filter")
    use_cache     = not getattr(args, "no_cache", False)
    no_plots      = getattr(args, "no_plots", False)
    no_genes      = getattr(args, "no_genes", False)
    job_timeout   = int(_get("job_timeout", 900))
    poll_interval = int(_get("poll_interval", 15))
    design_caveats = params.get("design_caveats", [])

    # ── Load experiment IDs ───────────────────────────────────────────────── #
    if args.group_file:
        group_col = _get("group_column", "group")
        experiments_a, experiments_b = load_two_group_file(
            args.group_file, group_column=group_col
        )
    elif args.group_a and args.group_b:
        experiments_a = load_ids_from_file(args.group_a)
        experiments_b = load_ids_from_file(args.group_b)
    elif "experiments_a" in params and "experiments_b" in params:
        experiments_a = params["experiments_a"]
        experiments_b = params["experiments_b"]
    else:
        print(
            "ERROR: Provide experiment IDs via --group-a/--group-b, "
            "--group-file, or a YAML config.",
            file=sys.stderr,
        )
        return 1

    # ── Run pipeline ──────────────────────────────────────────────────────── #
    try:
        run_diff_workflow(
            experiments_a=experiments_a,
            experiments_b=experiments_b,
            genome=genome,
            analysis_type=analysis_type,
            title=title,
            description_a=description_a,
            description_b=description_b,
            min_score=min_score,
            min_region_size=min_region_size,
            direction_filter=direction_filter,
            design_caveats=design_caveats,
            output_dir=output_dir,
            use_cache=use_cache,
            generate_plots=not no_plots,
            annotate_genes=not no_genes,
            job_timeout=job_timeout,
            poll_interval=poll_interval,
        )
    except (ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


def cmd_example(args: argparse.Namespace) -> int:
    """Run a built-in example dataset."""
    from src.utils import setup_logging
    from src.diff_analysis import run_diff_workflow

    setup_logging()

    dataset = args.dataset

    datasets = {
        "tp53_k562_dmso_vs_daunorubicin": {
            "experiments_a": [
                "SRX5865959", "SRX5865960", "SRX5865961", "SRX5865962",
                "SRX5865963", "SRX5865964", "SRX5865965", "SRX5865966",
            ],
            "experiments_b": [
                "SRX5865967", "SRX5865968", "SRX5865969", "SRX5865970",
                "SRX5865971", "SRX5865972", "SRX5865973", "SRX5865974",
            ],
            "genome":        "hg38",
            "analysis_type": "diffbind",
            "description_a": "K562_DMSO",
            "description_b": "K562_Daunorubicin",
            "title":         "TP53_K562_DMSO_vs_Daunorubicin",
            "design_caveats": [
                "Each group contains 8 different TP53 genotypes (WT, KO, and 6 hotspot "
                "mutants). edgeR treats all within-group experiments as replicates, so "
                "genotype-specific effects are averaged out."
            ],
        },
        "tp53_molm13_vs_k562": {
            "experiments_a": [
                "SRX5865975", "SRX5865976", "SRX5865977", "SRX5865982",
            ],
            "experiments_b": [
                "SRX5865959", "SRX5865960", "SRX5865961", "SRX5865966",
            ],
            "genome":        "hg38",
            "analysis_type": "diffbind",
            "description_a": "MOLM-13",
            "description_b": "K-562",
            "title":         "TP53_MOLM13_vs_K562",
            "design_caveats": [
                "Cross-cell-line comparison (MOLM-13 male vs K-562 female). "
                "Expected chrY sex-chromosome confound — demonstrates automated QC detection."
            ],
        },
    }

    if dataset not in datasets:
        print(
            f"ERROR: Unknown dataset '{dataset}'. "
            f"Available: {list(datasets.keys())}",
            file=sys.stderr,
        )
        return 1

    cfg = datasets[dataset]
    output_dir = args.output_dir or f"results/{dataset}"

    try:
        run_diff_workflow(
            experiments_a=cfg["experiments_a"],
            experiments_b=cfg["experiments_b"],
            genome=cfg["genome"],
            analysis_type=cfg["analysis_type"],
            description_a=cfg["description_a"],
            description_b=cfg["description_b"],
            title=cfg["title"],
            design_caveats=cfg.get("design_caveats", []),
            output_dir=output_dir,
        )
    except (ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


def cmd_cache(args: argparse.Namespace) -> int:
    """Handle the ``cache`` sub-command."""
    from src.cache import AnalysisCache
    cache = AnalysisCache()

    if args.cache_action == "list":
        keys = cache.list_keys()
        if keys:
            print(f"Cached analyses ({len(keys)}):")
            for k in keys:
                print(f"  {k}")
        else:
            print("Cache is empty.")
    elif args.cache_action == "clear":
        n = cache.clear()
        print(f"Deleted {n} cached analysis/analyses.")
    return 0


# ─────────────────────────── Parser construction ───────────────────────────── #


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="chipatlas-diff",
        description="ChIP-Atlas differential binding / methylation analysis pipeline.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--log-level", default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity (default: INFO)",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # ── run ────────────────────────────────────────────────────────────────── #
    p_run = sub.add_parser("run", help="Submit and run a differential analysis")
    p_run.add_argument("--config",    metavar="YAML", help="YAML config file")
    p_run.add_argument("--group-a",   metavar="FILE", help="Group A IDs (one per line)")
    p_run.add_argument("--group-b",   metavar="FILE", help="Group B IDs (one per line)")
    p_run.add_argument(
        "--group-file", metavar="FILE",
        help="Single file with two groups (requires --group-column)",
    )
    p_run.add_argument("--group-column", default="group", metavar="COL")
    p_run.add_argument("--genome",        default="hg38")
    p_run.add_argument("--analysis-type", default="diffbind",
                       choices=["diffbind", "dmr"])
    p_run.add_argument("--title",         default="chipatlas_diff")
    p_run.add_argument("--description-a", default="group_A")
    p_run.add_argument("--description-b", default="group_B")
    p_run.add_argument("--output-dir",    default="diff_analysis_results",
                       metavar="DIR")
    p_run.add_argument("--min-score",       type=int, default=0)
    p_run.add_argument("--min-region-size", type=int, default=0, metavar="BP")
    p_run.add_argument("--direction-filter",
                       choices=["A_enriched", "B_enriched"], default=None)
    p_run.add_argument("--no-cache",  action="store_true",
                       help="Force API re-submission even if cached")
    p_run.add_argument("--no-plots",  action="store_true",
                       help="Skip visualisation generation")
    p_run.add_argument("--no-genes",  action="store_true",
                       help="Skip UCSC gene annotation")
    p_run.add_argument("--job-timeout",   type=int, default=900, metavar="SEC")
    p_run.add_argument("--poll-interval", type=int, default=15,  metavar="SEC")
    p_run.set_defaults(func=cmd_run)

    # ── example ────────────────────────────────────────────────────────────── #
    p_ex = sub.add_parser("example", help="Run a built-in example dataset")
    p_ex.add_argument(
        "--dataset",
        default="tp53_k562_dmso_vs_daunorubicin",
        choices=["tp53_k562_dmso_vs_daunorubicin", "tp53_molm13_vs_k562"],
        help="Which built-in dataset to run",
    )
    p_ex.add_argument("--output-dir", default=None, metavar="DIR")
    p_ex.set_defaults(func=cmd_example)

    # ── cache ──────────────────────────────────────────────────────────────── #
    p_cache = sub.add_parser("cache", help="Manage the local analysis cache")
    p_cache.add_argument(
        "cache_action", choices=["list", "clear"],
        help="'list' all cached keys or 'clear' the entire cache",
    )
    p_cache.set_defaults(func=cmd_cache)

    return parser


def main(argv=None) -> int:
    parser = build_parser()
    args   = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
