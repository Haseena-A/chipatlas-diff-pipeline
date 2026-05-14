"""
Unit tests for chipatlas-diff-pipeline.

Run with:
    pytest tests/ -v
"""

from __future__ import annotations

import io
import os
import sys
import textwrap
import tempfile
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

# ─────────────────────────── validation ────────────────────────────────────── #

from src.validation import (
    get_standard_chroms,
    is_valid_experiment_id,
    validate_experiment_list,
    validate_genome,
    validate_analysis_type,
    run_preflight_checks,
)


class TestValidation:
    def test_valid_srx(self):
        assert is_valid_experiment_id("SRX12345")

    def test_valid_gsm(self):
        assert is_valid_experiment_id("GSM987654")

    def test_invalid_id(self):
        assert not is_valid_experiment_id("NOTVALID")
        assert not is_valid_experiment_id("SRA12345")

    def test_validate_genome_ok(self):
        validate_genome("hg38")    # no exception
        validate_genome("mm10")

    def test_validate_genome_bad(self):
        with pytest.raises(ValueError, match="Unsupported genome"):
            validate_genome("hg99")

    def test_validate_analysis_type(self):
        validate_analysis_type("diffbind")
        validate_analysis_type("dmr")
        with pytest.raises(ValueError):
            validate_analysis_type("unknown")

    def test_experiment_list_too_short(self):
        with pytest.raises(ValueError, match="at least 2"):
            validate_experiment_list(["SRX001"], "A")

    def test_preflight_ok(self):
        run_preflight_checks(["SRX001", "SRX002"], ["SRX003", "SRX004"], "hg38", "diffbind")

    def test_standard_chroms_hg38(self):
        chroms = get_standard_chroms("hg38")
        assert "chr1" in chroms
        assert "chrX" in chroms
        assert "chrM" in chroms
        assert "chrUn_random" not in chroms

    def test_standard_chroms_dm6(self):
        chroms = get_standard_chroms("dm6")
        assert "chr2L" in chroms
        assert "chrM" in chroms

    def test_standard_chroms_fallback(self):
        chroms = get_standard_chroms("UNKNOWN_GENOME")
        assert "chr1" in chroms   # falls back to hg38


# ─────────────────────────── parser ────────────────────────────────────────── #

from src.parser import parse_bed_results, summarize_regions, _mean_counts


def _make_bed_file(rows: list[list], tmp_dir: str) -> str:
    """Write a minimal tab-separated BED to a temp file."""
    path = os.path.join(tmp_dir, "wabi_result.bed")
    with open(path, "w") as fh:
        for row in rows:
            fh.write("\t".join(str(x) for x in row) + "\n")
    return path


class TestParser:
    BED_ROWS = [
        ["chr1", 1000, 2000, "50.0,60.0", "10.0,15.0", 2.5,  0.001, 0.01],
        ["chr2", 5000, 6000, "5.0,8.0",   "30.0,35.0", -1.8, 0.003, 0.02],
        ["chrX", 1000, 1050, "20.0",       "2.0",       3.0,  0.0001, 0.001],
    ]

    def test_parse_basic(self, tmp_path):
        bed = _make_bed_file(self.BED_ROWS, str(tmp_path))
        df = parse_bed_results(bed)
        assert len(df) == 3
        assert "logFC" in df.columns
        assert "qvalue" in df.columns
        assert "direction" in df.columns
        assert "significant" in df.columns
        assert "region_size" in df.columns

    def test_direction_assignment(self, tmp_path):
        bed = _make_bed_file(self.BED_ROWS, str(tmp_path))
        df = parse_bed_results(bed)
        a_row = df[df["chrom"] == "chr1"].iloc[0]
        b_row = df[df["chrom"] == "chr2"].iloc[0]
        assert a_row["direction"] == "A_enriched"
        assert b_row["direction"] == "B_enriched"

    def test_significant_flag(self, tmp_path):
        bed = _make_bed_file(self.BED_ROWS, str(tmp_path))
        df = parse_bed_results(bed)
        assert df[df["qvalue"] < 0.05]["significant"].all()

    def test_region_size(self, tmp_path):
        bed = _make_bed_file(self.BED_ROWS, str(tmp_path))
        df = parse_bed_results(bed)
        chr1_row = df[df["chrom"] == "chr1"].iloc[0]
        assert chr1_row["region_size"] == 1000

    def test_mean_counts(self):
        assert _mean_counts("10.0,20.0,30.0") == pytest.approx(20.0)
        assert _mean_counts("bad_data") == 0.0
        assert _mean_counts("50.0") == 50.0

    def test_empty_bed_raises(self, tmp_path):
        path = str(tmp_path / "empty.bed")
        with open(path, "w") as fh:
            fh.write("")
        with pytest.raises(RuntimeError):
            parse_bed_results(path)

    def test_sorted_by_qvalue(self, tmp_path):
        bed = _make_bed_file(self.BED_ROWS, str(tmp_path))
        df = parse_bed_results(bed)
        assert df["qvalue"].is_monotonic_increasing


# ─────────────────────────── filters ───────────────────────────────────────── #

from src.filters import filter_standard_chromosomes, filter_min_region_size, filter_regions


def _sample_df() -> pd.DataFrame:
    return pd.DataFrame({
        "chrom":       ["chr1", "chr2", "chrUn_gl000220", "chrM", "chr3"],
        "chromStart":  [100,    200,    300,              400,    500],
        "chromEnd":    [600,    700,    800,              900,    520],
        "logFC":       [1.5,   -2.0,    0.5,              3.0,   -0.5],
        "qvalue":      [0.01,   0.001,  0.04,             0.001,  0.06],
        "direction":   ["A_enriched","B_enriched","A_enriched","A_enriched","B_enriched"],
        "significant": [True,  True,   True,             True,   False],
        "region_size": [500,   500,    500,              500,    20],
        "score":       [300,   400,    200,              500,    100],
    })


class TestFilters:
    def test_standard_chroms_removes_nonstandard(self):
        df = _sample_df()
        out = filter_standard_chromosomes(df, genome="hg38")
        assert "chrUn_gl000220" not in out["chrom"].values
        assert len(out) == 4

    def test_min_region_size(self):
        df = _sample_df()
        out = filter_min_region_size(df, min_size=100)
        assert all(out["region_size"] >= 100)

    def test_filter_regions_max_qvalue(self):
        df = _sample_df()
        out = filter_regions(df, max_qvalue=0.05)
        assert (out["qvalue"] < 0.05).all()

    def test_filter_regions_direction(self):
        df = _sample_df()
        out = filter_regions(df, direction_filter="A_enriched")
        assert (out["direction"] == "A_enriched").all()

    def test_filter_regions_min_logfc(self):
        df = _sample_df()
        out = filter_regions(df, min_logfc=1.0)
        assert (out["logFC"].abs() >= 1.0).all()

    def test_filter_regions_chromosomes(self):
        df = _sample_df()
        out = filter_regions(df, chromosomes=["chr1", "chr2"])
        assert set(out["chrom"]) == {"chr1", "chr2"}


# ─────────────────────────── qc_checks ─────────────────────────────────────── #

from src.qc_checks import run_qc_checks, _check_sample_size, _check_genomic_clusters


class TestQCChecks:
    def test_no_warnings_on_clean_data(self):
        df = pd.DataFrame({
            "chrom":       ["chr1"] * 5,
            "chromStart":  [100_000 * i for i in range(5)],
            "chromEnd":    [100_000 * i + 500 for i in range(5)],
            "logFC":       [1.5, -1.5, 2.0, -2.0, 1.0],
            "qvalue":      [0.001] * 5,
            "direction":   ["A_enriched", "B_enriched"] * 2 + ["A_enriched"],
            "significant": [True] * 5,
            "region_size": [500] * 5,
            "mean_count_a": [100.0] * 5,
            "mean_count_b": [50.0] * 5,
        })
        warnings = run_qc_checks(df, n_samples_a=5, n_samples_b=5)
        severity = [w["severity"] for w in warnings]
        assert "major" not in severity

    def test_low_sample_size_major(self):
        warnings = _check_sample_size(2, 2)
        assert len(warnings) == 1
        assert warnings[0]["severity"] == "major"

    def test_low_sample_size_minor(self):
        warnings = _check_sample_size(4, 4)
        assert len(warnings) == 1
        assert warnings[0]["severity"] == "minor"

    def test_adequate_sample_size(self):
        assert _check_sample_size(6, 6) == []

    def test_cluster_detection(self):
        # 12 significant regions within 1 Mb → should be flagged
        df = pd.DataFrame({
            "chrom":       ["chr1"] * 12,
            "chromStart":  list(range(0, 12 * 50_000, 50_000)),
            "chromEnd":    list(range(500, 12 * 50_000 + 500, 50_000)),
            "direction":   ["A_enriched"] * 12,
            "significant": [True] * 12,
            "mean_count_a": [100.0] * 12,
            "mean_count_b": [1.0] * 12,
        })
        warnings = _check_genomic_clusters(df)
        assert len(warnings) >= 1
        assert "genomic cluster" in warnings[0]["issue"]


# ─────────────────────────── cache ─────────────────────────────────────────── #

from src.cache import AnalysisCache, make_cache_key


class TestCache:
    def test_set_get(self, tmp_path):
        cache = AnalysisCache(cache_dir=str(tmp_path))
        key = "testkey"
        data = {"result": [1, 2, 3]}
        cache.set(key, data)
        loaded = cache.get(key)
        assert loaded == data

    def test_miss_returns_none(self, tmp_path):
        cache = AnalysisCache(cache_dir=str(tmp_path))
        assert cache.get("nonexistent") is None

    def test_delete(self, tmp_path):
        cache = AnalysisCache(cache_dir=str(tmp_path))
        cache.set("k", {"x": 1})
        assert cache.delete("k") is True
        assert cache.get("k") is None

    def test_list_keys(self, tmp_path):
        cache = AnalysisCache(cache_dir=str(tmp_path))
        cache.set("a", {})
        cache.set("b", {})
        assert set(cache.list_keys()) == {"a", "b"}

    def test_cache_key_deterministic(self):
        k1 = make_cache_key(["SRX001", "SRX002"], ["SRX003"], "hg38", "diffbind")
        k2 = make_cache_key(["SRX002", "SRX001"], ["SRX003"], "hg38", "diffbind")
        assert k1 == k2   # order-independent

    def test_cache_key_different_genomes(self):
        k1 = make_cache_key(["SRX001", "SRX002"], ["SRX003"], "hg38", "diffbind")
        k2 = make_cache_key(["SRX001", "SRX002"], ["SRX003"], "mm10", "diffbind")
        assert k1 != k2


# ─────────────────────────── utils ─────────────────────────────────────────── #

from src.utils import load_ids_from_file, load_two_group_file


class TestUtils:
    def test_load_plain_text(self, tmp_path):
        f = tmp_path / "ids.txt"
        f.write_text("SRX001\n# comment\nSRX002\n")
        ids = load_ids_from_file(str(f))
        assert ids == ["SRX001", "SRX002"]

    def test_load_csv(self, tmp_path):
        f = tmp_path / "ids.csv"
        f.write_text("experiment_id,extra\nSRX001,foo\nSRX002,bar\n")
        ids = load_ids_from_file(str(f))
        assert "SRX001" in ids
        assert "SRX002" in ids

    def test_load_two_group_file(self, tmp_path):
        f = tmp_path / "meta.csv"
        f.write_text(
            "experiment_id,group\n"
            "SRX001,control\n"
            "SRX002,control\n"
            "SRX003,treated\n"
            "SRX004,treated\n"
        )
        a, b = load_two_group_file(str(f))
        assert sorted(a) == ["SRX001", "SRX002"]
        assert sorted(b) == ["SRX003", "SRX004"]

    def test_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            load_ids_from_file("/nonexistent/path.txt")


# ─────────────────────────── api_client (mocked) ──────────────────────────── #

from src.api_client import validate_inputs, validate_experiment_ids


class TestAPIClientValidation:
    def test_validate_inputs_ok(self):
        validate_inputs(["SRX001", "SRX002"], ["SRX003", "SRX004"], "hg38", "diffbind")

    def test_validate_inputs_bad_genome(self):
        with pytest.raises(ValueError, match="Supported assemblies"):
            validate_inputs(["SRX001", "SRX002"], ["SRX003", "SRX004"], "hg99", "diffbind")

    def test_validate_inputs_too_few(self):
        with pytest.raises(ValueError, match="at least 2"):
            validate_inputs(["SRX001"], ["SRX002", "SRX003"], "hg38", "diffbind")
