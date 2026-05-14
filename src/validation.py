"""
Input validation for ChIP-Atlas differential analysis.

Provides genome-aware chromosome sets, experiment ID validation,
and request pre-flight checks.
"""

from __future__ import annotations

import re
from typing import Dict, List, Optional, Set

# ─────────────────────────── Genome chromosome sets ───────────────────────── #

GENOME_STANDARD_CHROMS: Dict[str, Set[str]] = {
    # Human
    "hg38": {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y", "M"]},
    "hg19": {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y", "M"]},
    # Mouse
    "mm10": {f"chr{c}" for c in list(range(1, 20)) + ["X", "Y", "M"]},
    "mm9":  {f"chr{c}" for c in list(range(1, 20)) + ["X", "Y", "M"]},
    # Rat
    "rn6":  {f"chr{c}" for c in list(range(1, 21)) + ["X", "Y", "M"]},
    # Drosophila
    "dm6": {"chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX", "chrY", "chrM"},
    "dm3": {"chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX", "chrY", "chrM"},
    # C. elegans
    "ce11": {"chrI", "chrII", "chrIII", "chrIV", "chrV", "chrX", "chrM"},
    "ce10": {"chrI", "chrII", "chrIII", "chrV", "chrV", "chrX", "chrM"},
    # Yeast
    "sacCer3": {
        f"chr{r}" for r in [
            "I", "II", "III", "IV", "V", "VI", "VII", "VIII",
            "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI", "M",
        ]
    },
}

EXPERIMENT_ID_RE = re.compile(r"^(SRX|ERX|DRX|GSM)\d+$", re.IGNORECASE)

VALID_GENOMES = list(GENOME_STANDARD_CHROMS.keys())
VALID_ANALYSIS_TYPES = ["diffbind", "dmr"]


# ─────────────────────────── Public helpers ────────────────────────────────── #


def get_standard_chroms(genome: Optional[str] = None) -> Set[str]:
    """
    Return the standard chromosome set for a genome assembly.

    Falls back to human (hg38) if ``genome`` is None or unrecognised.

    Parameters
    ----------
    genome : str or None
        Assembly identifier (e.g. 'hg38', 'mm10').

    Returns
    -------
    set of str
        Standard chromosome names (e.g. ``{'chr1', ..., 'chrX', 'chrM'}``).
    """
    if genome and genome in GENOME_STANDARD_CHROMS:
        return GENOME_STANDARD_CHROMS[genome]
    return GENOME_STANDARD_CHROMS["hg38"]


def is_valid_experiment_id(exp_id: str) -> bool:
    """Return True if *exp_id* matches the SRX/ERX/DRX/GSM pattern."""
    return bool(EXPERIMENT_ID_RE.match(exp_id))


def validate_experiment_list(ids: List[str], group_label: str = "") -> None:
    """
    Raise ValueError if the experiment list is too short; warn on bad IDs.

    Parameters
    ----------
    ids : list of str
        Experiment accessions.
    group_label : str
        Label used in error messages (e.g. 'A' or 'B').

    Raises
    ------
    ValueError
        If fewer than 2 IDs are provided.
    """
    if len(ids) < 2:
        raise ValueError(
            f"Group {group_label}: at least 2 experiment IDs are required "
            f"(got {len(ids)})."
        )
    invalid = [i for i in ids if not is_valid_experiment_id(i)]
    if invalid:
        print(
            f"   Warning – Group {group_label}: {len(invalid)} ID(s) don't "
            f"match SRX/ERX/DRX/GSM format: {invalid[:5]}"
        )


def validate_genome(genome: str) -> None:
    """
    Raise ValueError if *genome* is not a supported assembly.

    Parameters
    ----------
    genome : str
        Assembly name to check.

    Raises
    ------
    ValueError
        If the assembly is not supported.
    """
    if genome not in VALID_GENOMES:
        raise ValueError(
            f"Unsupported genome '{genome}'. "
            f"Supported assemblies: {VALID_GENOMES}"
        )


def validate_analysis_type(analysis_type: str) -> None:
    """
    Raise ValueError if *analysis_type* is not supported.

    Parameters
    ----------
    analysis_type : str
        Analysis type string.

    Raises
    ------
    ValueError
    """
    if analysis_type not in VALID_ANALYSIS_TYPES:
        raise ValueError(
            f"Unsupported analysis_type '{analysis_type}'. "
            f"Supported: {VALID_ANALYSIS_TYPES}"
        )


def run_preflight_checks(
    experiments_a: List[str],
    experiments_b: List[str],
    genome: str,
    analysis_type: str,
) -> None:
    """
    Run all pre-submission validation checks in one call.

    Raises
    ------
    ValueError
        On any validation failure.
    """
    validate_experiment_list(experiments_a, "A")
    validate_experiment_list(experiments_b, "B")
    validate_genome(genome)
    validate_analysis_type(analysis_type)
