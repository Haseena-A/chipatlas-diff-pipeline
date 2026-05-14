"""
ChIP-Atlas Differential Analysis API Client.

Handles job submission, status polling, and result retrieval from the
official ChIP-Atlas WABI API:

    POST https://dtn1.ddbj.nig.ac.jp/wabi/chipatlas/

Reference: Zou et al. 2024 (NAR), https://chip-atlas.org
"""

from __future__ import annotations

import io
import logging
import os
import re
import time
import zipfile
from typing import Dict, List, Optional

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

logger = logging.getLogger(__name__)

# ─────────────────────────── Constants ────────────────────────────────────── #

WABI_API_URL = "https://dtn1.ddbj.nig.ac.jp/wabi/chipatlas/"

VALID_ANALYSIS_TYPES: Dict[str, str] = {
    "diffbind": "Differential Peak Regions (ChIP/ATAC/DNase-seq via edgeR)",
    "dmr": "Differentially Methylated Regions (Bisulfite-seq via metilene)",
}

VALID_GENOMES: List[str] = [
    "hg38", "hg19",
    "mm10", "mm9",
    "rn6",
    "dm6", "dm3",
    "ce11", "ce10",
    "sacCer3",
]

EXPERIMENT_ID_PATTERN = re.compile(r"^(SRX|ERX|DRX|GSM)\d+$", re.IGNORECASE)

# ─────────────────────────── Session helpers ───────────────────────────────── #


def _build_session(retries: int = 3, backoff: float = 1.0) -> requests.Session:
    """Return a requests.Session with retry logic."""
    session = requests.Session()
    retry = Retry(
        total=retries,
        backoff_factor=backoff,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET", "POST"],
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session


# ─────────────────────────── Validation ───────────────────────────────────── #


def validate_experiment_ids(ids: List[str], group_name: str = "") -> List[str]:
    """
    Warn about experiment IDs that don't match expected format.

    Parameters
    ----------
    ids : list of str
        Accession IDs to validate.
    group_name : str
        Label for log messages (e.g. 'A' or 'B').

    Returns
    -------
    list of str
        The original list (invalid IDs are warned, not dropped).
    """
    invalid = [i for i in ids if not EXPERIMENT_ID_PATTERN.match(i)]
    if invalid:
        logger.warning(
            "Group %s: %d ID(s) don't match SRX/ERX/DRX/GSM pattern: %s",
            group_name,
            len(invalid),
            invalid[:5],
        )
    return ids


def validate_inputs(
    experiments_a: List[str],
    experiments_b: List[str],
    genome: str,
    analysis_type: str,
) -> None:
    """
    Validate all parameters before API submission.

    Raises
    ------
    ValueError
        On any invalid parameter.
    """
    if len(experiments_a) < 2:
        raise ValueError("Group A requires at least 2 experiment IDs.")
    if len(experiments_b) < 2:
        raise ValueError("Group B requires at least 2 experiment IDs.")
    if genome not in VALID_GENOMES:
        raise ValueError(
            f"Invalid genome '{genome}'. Supported assemblies: {VALID_GENOMES}"
        )
    if analysis_type not in VALID_ANALYSIS_TYPES:
        raise ValueError(
            f"Invalid analysis_type '{analysis_type}'. "
            f"Supported: {list(VALID_ANALYSIS_TYPES.keys())}"
        )
    validate_experiment_ids(experiments_a, "A")
    validate_experiment_ids(experiments_b, "B")


# ─────────────────────────── API calls ────────────────────────────────────── #


def submit_diff_job(
    experiments_a: List[str],
    experiments_b: List[str],
    genome: str = "hg38",
    analysis_type: str = "diffbind",
    title: str = "diff_analysis",
    description_a: str = "group_A",
    description_b: str = "group_B",
    timeout: int = 60,
) -> str:
    """
    Submit a differential analysis job to the ChIP-Atlas WABI API.

    Parameters
    ----------
    experiments_a : list of str
        SRX/ERX/DRX/GSM accessions for group A (≥ 2 required).
    experiments_b : list of str
        SRX/ERX/DRX/GSM accessions for group B (≥ 2 required).
    genome : str
        UCSC genome assembly (e.g. 'hg38', 'mm10').
    analysis_type : str
        'diffbind' (edgeR DPR) or 'dmr' (metilene DMR).
    title : str
        Short analysis title stored in job metadata.
    description_a : str
        Human-readable label for group A.
    description_b : str
        Human-readable label for group B.
    timeout : int
        HTTP request timeout in seconds.

    Returns
    -------
    str
        ChIP-Atlas request ID for polling.

    Raises
    ------
    ValueError
        If inputs are invalid.
    RuntimeError
        If submission fails or response cannot be parsed.
    """
    validate_inputs(experiments_a, experiments_b, genome, analysis_type)

    payload = {
        "address": "",
        "format": "text",
        "result": "www",
        "genome": genome,
        "antigenClass": analysis_type,
        "cellClass": "empty",
        "threshold": "50",
        "typeA": "srx",
        "bedAFile": "\n".join(experiments_a),
        "typeB": "srx",
        "bedBFile": "\n".join(experiments_b),
        "permTime": "1",
        "title": title,
        "descriptionA": description_a,
        "descriptionB": description_b,
    }

    logger.info(
        "Submitting %s job to ChIP-Atlas (genome=%s, n_A=%d, n_B=%d) …",
        analysis_type, genome, len(experiments_a), len(experiments_b),
    )
    session = _build_session()

    try:
        resp = session.post(WABI_API_URL, data=payload, timeout=timeout)
        resp.raise_for_status()
    except requests.exceptions.RequestException as exc:
        raise RuntimeError(
            f"ChIP-Atlas API unreachable: {exc}. "
            "Check internet connection and https://chip-atlas.org status."
        ) from exc

    # Try JSON first
    try:
        body = resp.json()
        if "requestId" in body:
            rid = body["requestId"]
            logger.info("Job submitted. Request ID: %s", rid)
            return rid
        if err := body.get("error-message", body.get("Message", "")):
            raise RuntimeError(f"API returned an error: {err}")
    except (ValueError, KeyError):
        pass

    # Plain-text line parsing
    for line in resp.text.strip().splitlines():
        if line.lower().startswith("requestid"):
            parts = re.split(r"[:\t]\s*", line, maxsplit=1)
            if len(parts) == 2 and parts[1].strip():
                rid = parts[1].strip()
                logger.info("Job submitted. Request ID: %s", rid)
                return rid

    raise RuntimeError(
        f"Could not parse a request ID from the API response:\n{resp.text[:500]}"
    )


def poll_job_status(
    request_id: str,
    timeout: int = 900,
    poll_interval: int = 15,
) -> str:
    """
    Block until the ChIP-Atlas job finishes or times out.

    Parameters
    ----------
    request_id : str
        ID returned by :func:`submit_diff_job`.
    timeout : int
        Maximum wall-clock seconds to wait (default 900 = 15 min).
    poll_interval : int
        Seconds between status checks.

    Returns
    -------
    str
        Final status string ('finished').

    Raises
    ------
    RuntimeError
        If the job errors or the timeout is exceeded.
    """
    url = f"{WABI_API_URL}{request_id}?info=status"
    session = _build_session()
    t0 = time.monotonic()

    while True:
        elapsed = time.monotonic() - t0
        if elapsed > timeout:
            raise RuntimeError(
                f"Job {request_id} timed out after {timeout}s. "
                "Try increasing --timeout."
            )

        try:
            resp = session.get(url, timeout=30)
            resp.raise_for_status()
            status: Optional[str] = None
            for line in resp.text.strip().splitlines():
                if line.lower().startswith("status"):
                    parts = re.split(r"[:\t]\s*", line, maxsplit=1)
                    status = parts[1].strip() if len(parts) == 2 else None
                    break

            if status == "finished":
                logger.info("Job %s finished after %.0fs.", request_id, elapsed)
                return "finished"
            if status in ("error", "failed"):
                raise RuntimeError(
                    f"ChIP-Atlas job {request_id} failed. "
                    f"Response: {resp.text[:300]}"
                )

            mins = elapsed / 60
            logger.info("  [%.1f min] Job running (status=%s) …", mins, status)
            print(f"   Job running … ({mins:.1f} min elapsed, status={status})")

        except requests.exceptions.RequestException as exc:
            logger.warning("Poll request failed (%s), retrying …", exc)
            print(f"   Warning: poll failed ({exc}), retrying …")

        time.sleep(poll_interval)


def retrieve_zip_results(
    request_id: str,
    output_dir: str,
    http_timeout: int = 300,
) -> Dict[str, str]:
    """
    Download and extract the ZIP result archive from ChIP-Atlas.

    Parameters
    ----------
    request_id : str
        Completed job ID.
    output_dir : str
        Parent directory; extracted files go into ``<output_dir>/raw_results/``.
    http_timeout : int
        Download timeout in seconds.

    Returns
    -------
    dict
        ``{'bed': path, 'igv_bed': path, 'log': path, 'igv_xml': path}``
        (absent keys means not found in ZIP).

    Raises
    ------
    RuntimeError
        If download, ZIP extraction, or BED lookup fails.
    """
    url = f"{WABI_API_URL}{request_id}?info=result&format=zip"
    raw_dir = os.path.join(output_dir, "raw_results")
    os.makedirs(raw_dir, exist_ok=True)

    logger.info("Downloading results ZIP from %s …", url)
    session = _build_session(retries=2)

    try:
        resp = session.get(url, timeout=http_timeout)
        resp.raise_for_status()
        if not resp.content:
            raise RuntimeError("API returned an empty response body.")
    except requests.exceptions.RequestException as exc:
        raise RuntimeError(f"Failed to download results: {exc}") from exc

    try:
        with zipfile.ZipFile(io.BytesIO(resp.content)) as zf:
            zf.extractall(raw_dir)
    except zipfile.BadZipFile as exc:
        raise RuntimeError(
            f"API returned invalid ZIP data. Job may still be processing. "
            f"Request ID: {request_id}"
        ) from exc

    extracted: Dict[str, str] = {}
    for root, _dirs, files in os.walk(raw_dir):
        for fname in files:
            fpath = os.path.join(root, fname)
            if fname.endswith(".igv.bed"):
                extracted["igv_bed"] = fpath
            elif fname.endswith(".bed") and "igv" not in fname:
                extracted["bed"] = fpath
            elif fname.endswith(".log"):
                extracted["log"] = fpath
            elif fname.endswith(".xml"):
                extracted["igv_xml"] = fpath

    if "bed" not in extracted:
        all_files = [f for _, _, fs in os.walk(raw_dir) for f in fs]
        raise RuntimeError(
            f"No BED file found in ZIP archive. Extracted: {all_files}"
        )

    logger.info("Extracted %d file(s) to %s", len(extracted), raw_dir)
    return extracted
