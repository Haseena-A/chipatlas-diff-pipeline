#!/usr/bin/env python3
"""
Example: TP53 ChIP-seq differential binding in K562 cells.

Compares DMSO (vehicle) vs Daunorubicin (chemotherapy) to identify
genomic regions where DNA damage alters TP53 occupancy.

Dataset: GSE131484 (Bhatt et al.)
Genome:  hg38
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from src.diff_analysis import run_diff_workflow

# ── Experiment groups ─────────────────────────────────────────────────────── #

GROUP_A = [               # K562 cells + DMSO (vehicle control)
    "SRX5865959",         # R175H mutant
    "SRX5865960",         # Y220C mutant
    "SRX5865961",         # M237I mutant
    "SRX5865962",         # R248Q mutant
    "SRX5865963",         # R273H mutant
    "SRX5865964",         # R282W mutant
    "SRX5865965",         # KO
    "SRX5865966",         # WT
]

GROUP_B = [               # K562 cells + Daunorubicin (100 nM, 24 h)
    "SRX5865967",
    "SRX5865968",
    "SRX5865969",
    "SRX5865970",
    "SRX5865971",
    "SRX5865972",
    "SRX5865973",
    "SRX5865974",
]

DESIGN_CAVEATS = [
    "Each group contains 8 different TP53 genotypes (WT, KO, and 6 hotspot "
    "mutants). edgeR treats all within-group experiments as replicates, so "
    "genotype-specific effects are averaged out. The analysis captures the "
    "average DMSO-vs-Daunorubicin effect across genotypes."
]

# ── Run pipeline ──────────────────────────────────────────────────────────── #

if __name__ == "__main__":
    results = run_diff_workflow(
        experiments_a  = GROUP_A,
        experiments_b  = GROUP_B,
        genome         = "hg38",
        analysis_type  = "diffbind",
        title          = "TP53_K562_DMSO_vs_Daunorubicin",
        description_a  = "K562_DMSO",
        description_b  = "K562_Daunorubicin",
        design_caveats = DESIGN_CAVEATS,
        output_dir     = "results/tp53_k562_dmso_vs_daunorubicin",
        use_cache      = True,
        generate_plots = True,
        annotate_genes = True,
    )

    # ── Quick downstream access ───────────────────────────────────────────── #
    df  = results["diff_regions"]
    sig = df[df["significant"]]

    print(f"\nDownstream summary:")
    print(f"  Significant regions: {len(sig):,}")
    print(f"  K562_DMSO-enriched:         {int((sig['direction']=='A_enriched').sum()):,}")
    print(f"  K562_Daunorubicin-enriched: {int((sig['direction']=='B_enriched').sum()):,}")

    if len(sig):
        top = sig.head(5)[["chrom", "chromStart", "chromEnd", "logFC", "qvalue"]]
        print(f"\nTop 5 significant regions:\n{top.to_string(index=False)}")
