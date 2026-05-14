"""
Utility helpers for the ChIP-Atlas differential analysis pipeline.

Includes logging setup, output directory management, and experiment ID
loading from plain-text / CSV / TSV files.
"""

from __future__ import annotations

import logging
import os
import re
import sys
from typing import List, Optional, Tuple

import pandas as pd

EXPERIMENT_ID_RE = re.compile(r"^(SRX|ERX|DRX|GSM)\d+$", re.IGNORECASE)


# ─────────────────────────── Logging ───────────────────────────────────────── #


def setup_logging(level: str = "INFO", log_file: Optional[str] = None) -> None:
    """
    Configure root logger with a consistent format.

    Parameters
    ----------
    level : str
        Logging level (``'DEBUG'``, ``'INFO'``, ``'WARNING'``, etc.).
    log_file : str or None
        If provided, logs are also written to this file path.
    """
    fmt = "%(asctime)s  %(levelname)-8s  %(name)s  %(message)s"
    handlers: List[logging.Handler] = [logging.StreamHandler(sys.stdout)]
    if log_file:
        os.makedirs(os.path.dirname(log_file) or ".", exist_ok=True)
        handlers.append(logging.FileHandler(log_file))

    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format=fmt,
        handlers=handlers,
        force=True,
    )


# ─────────────────────────── Directory management ─────────────────────────── #


def ensure_output_dir(path: str) -> str:
    """
    Create *path* (and all parents) if it does not exist.

    Returns
    -------
    str
        The resolved absolute path.
    """
    path = os.path.abspath(path)
    os.makedirs(path, exist_ok=True)
    return path


# ─────────────────────────── Experiment ID loading ─────────────────────────── #


def load_ids_from_file(
    file_path: str,
    column_name: Optional[str] = None,
) -> List[str]:
    """
    Load experiment IDs from a plain-text, CSV, or TSV file.

    Supported formats
    -----------------
    * **Plain text** – one accession per line (lines starting with ``#`` ignored)
    * **CSV / TSV**  – auto-detects or uses *column_name*

    Parameters
    ----------
    file_path : str
        Path to the input file.
    column_name : str or None
        Column name for the experiment IDs (auto-detected if None).

    Returns
    -------
    list of str
        Experiment accessions (may include invalid-format IDs with a warning).

    Raises
    ------
    FileNotFoundError
        If *file_path* does not exist.
    ValueError
        If no IDs can be found in the file.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")

    ext = os.path.splitext(file_path)[1].lower()
    if ext == ".txt":
        ids = _load_plain_text(file_path)
    elif ext in (".csv", ".tsv", ".tab"):
        sep = "," if ext == ".csv" else "\t"
        ids = _load_delimited(file_path, sep=sep, column_name=column_name)
    else:
        # Guess
        try:
            ids = _load_delimited(file_path, sep=",", column_name=column_name)
        except Exception:
            try:
                ids = _load_delimited(file_path, sep="\t", column_name=column_name)
            except Exception:
                ids = _load_plain_text(file_path)

    invalid = [i for i in ids if not EXPERIMENT_ID_RE.match(i)]
    if invalid:
        print(
            f"   Warning: {len(invalid)} ID(s) don't match SRX/ERX/DRX/GSM format: "
            f"{invalid[:5]}"
        )
    if not ids:
        raise ValueError(f"No experiment IDs found in {file_path}")
    print(f"✓ Loaded {len(ids)} experiment ID(s) from {file_path}")
    return ids


def load_two_group_file(
    file_path: str,
    group_column: str = "group",
    id_column: Optional[str] = None,
) -> Tuple[List[str], List[str]]:
    """
    Load two experiment groups from a single CSV/TSV with a group column.

    The file must have exactly two unique values in *group_column*.
    Groups are assigned alphabetically: first → A, second → B.

    Parameters
    ----------
    file_path : str
        Path to the two-group metadata file.
    group_column : str
        Column name containing group labels (default ``'group'``).
    id_column : str or None
        Column name for experiment IDs (auto-detected if None).

    Returns
    -------
    tuple of (list, list)
        ``(experiments_a, experiments_b)``

    Raises
    ------
    ValueError
        If the file doesn't have exactly 2 groups.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")

    ext = os.path.splitext(file_path)[1].lower()
    sep = "," if ext == ".csv" else "\t"
    df  = pd.read_csv(file_path, sep=sep)

    if group_column not in df.columns:
        raise ValueError(
            f"Group column '{group_column}' not found. "
            f"Available: {list(df.columns)}"
        )

    if id_column is None:
        _candidates = [
            "experiment_id", "srx", "srx_id", "SRX",
            "accession", "sample_id", "id", "ID",
        ]
        for c in _candidates:
            if c in df.columns:
                id_column = c
                break
        if id_column is None:
            id_column = next(c for c in df.columns if c != group_column)

    groups = sorted(df[group_column].unique())
    if len(groups) != 2:
        raise ValueError(
            f"Expected exactly 2 groups, found {len(groups)}: {groups}"
        )
    g_a, g_b = groups
    ids_a = df[df[group_column] == g_a][id_column].astype(str).tolist()
    ids_b = df[df[group_column] == g_b][id_column].astype(str).tolist()
    print(f"✓ Loaded groups: '{g_a}' ({len(ids_a)} IDs) → A, '{g_b}' ({len(ids_b)} IDs) → B")
    return ids_a, ids_b


# ─────────────────────────── Private helpers ───────────────────────────────── #


def _load_plain_text(file_path: str) -> List[str]:
    ids: List[str] = []
    with open(file_path) as fh:
        for line in fh:
            s = line.strip()
            if s and not s.startswith("#"):
                ids.append(s)
    if not ids:
        raise ValueError(f"No experiment IDs found in {file_path}")
    return ids


def _load_delimited(
    file_path: str,
    sep: str = ",",
    column_name: Optional[str] = None,
) -> List[str]:
    df = pd.read_csv(file_path, sep=sep)
    if column_name is None:
        _candidates = [
            "experiment_id", "srx", "srx_id", "SRX",
            "accession", "sample_id", "id", "ID",
        ]
        for c in _candidates:
            if c in df.columns:
                column_name = c
                break
        if column_name is None:
            column_name = df.columns[0]
            print(f"   Auto-detected ID column: '{column_name}'")

    if column_name not in df.columns:
        raise ValueError(
            f"Column '{column_name}' not found. Available: {list(df.columns)}"
        )

    ids = df[column_name].astype(str).tolist()
    ids = [x for x in ids if x not in ("nan", "")]
    if not ids:
        raise ValueError(f"No IDs found in column '{column_name}'")
    return ids
