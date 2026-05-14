"""
Simple disk-backed cache for completed ChIP-Atlas analyses.

Stores results keyed on a deterministic hash of the request parameters.
Prevents duplicate API submissions for identical inputs.

Cache directory: ``~/.chipatlas_cache/`` (configurable).
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import pickle
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)

DEFAULT_CACHE_DIR = os.path.expanduser("~/.chipatlas_cache")


# ─────────────────────────── Public API ────────────────────────────────────── #


class AnalysisCache:
    """
    Key-value cache backed by pickle files on disk.

    Parameters
    ----------
    cache_dir : str
        Directory where cache entries are stored.
    """

    def __init__(self, cache_dir: str = DEFAULT_CACHE_DIR) -> None:
        self.cache_dir = cache_dir
        os.makedirs(cache_dir, exist_ok=True)
        logger.debug("Cache directory: %s", cache_dir)

    # ── Core operations ─────────────────────────────────────────────────── #

    def get(self, key: str) -> Optional[Dict]:
        """
        Retrieve a cached result.

        Parameters
        ----------
        key : str
            Cache key (e.g. from :func:`make_cache_key`).

        Returns
        -------
        dict or None
            Cached results dict, or None if not found.
        """
        path = self._path(key)
        if not os.path.exists(path):
            return None
        try:
            with open(path, "rb") as fh:
                result = pickle.load(fh)
            logger.info("Cache hit: %s", key)
            return result
        except (pickle.UnpicklingError, EOFError, Exception) as exc:
            logger.warning("Cache read failed for %s: %s – ignoring.", key, exc)
            return None

    def set(self, key: str, value: Dict) -> None:
        """
        Store a result in the cache.

        Parameters
        ----------
        key : str
            Cache key.
        value : dict
            Results dict to cache.
        """
        path = self._path(key)
        try:
            with open(path, "wb") as fh:
                pickle.dump(value, fh)
            logger.info("Cached result under key: %s", key)
        except Exception as exc:
            logger.warning("Cache write failed for %s: %s", key, exc)

    def delete(self, key: str) -> bool:
        """
        Remove a cached entry.

        Returns
        -------
        bool
            True if the entry existed and was removed.
        """
        path = self._path(key)
        if os.path.exists(path):
            os.remove(path)
            logger.info("Deleted cache entry: %s", key)
            return True
        return False

    def list_keys(self) -> List[str]:
        """Return all cache keys currently on disk."""
        return [
            f[:-4]  # strip .pkl
            for f in os.listdir(self.cache_dir)
            if f.endswith(".pkl")
        ]

    def clear(self) -> int:
        """
        Delete all cached entries.

        Returns
        -------
        int
            Number of entries deleted.
        """
        count = 0
        for key in self.list_keys():
            self.delete(key)
            count += 1
        return count

    # ── Helpers ─────────────────────────────────────────────────────────── #

    def _path(self, key: str) -> str:
        return os.path.join(self.cache_dir, f"{key}.pkl")


def make_cache_key(
    experiments_a: List[str],
    experiments_b: List[str],
    genome: str,
    analysis_type: str,
) -> str:
    """
    Create a deterministic cache key from analysis parameters.

    The key is a 16-character hex digest of the sorted experiment IDs,
    genome, and analysis type — independent of list order.

    Parameters
    ----------
    experiments_a : list of str
    experiments_b : list of str
    genome : str
    analysis_type : str

    Returns
    -------
    str
        16-character hex string suitable for use as a filename stem.
    """
    payload = json.dumps({
        "a": sorted(experiments_a),
        "b": sorted(experiments_b),
        "genome": genome,
        "type": analysis_type,
    }, sort_keys=True)
    return hashlib.sha256(payload.encode()).hexdigest()[:16]
