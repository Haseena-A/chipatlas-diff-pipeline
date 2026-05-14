# chipatlas-diff-pipeline

[![Python](https://img.shields.io/badge/python-3.9%20%7C%203.10%20%7C%203.11%20%7C%203.12-blue?logo=python)](https://python.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![ChIP-Atlas](https://img.shields.io/badge/API-ChIP--Atlas-green)](https://chip-atlas.org)
[![edgeR](https://img.shields.io/badge/stats-edgeR%20DPR-orange)](https://bioconductor.org/packages/edgeR/)
[![Code style: ruff](https://img.shields.io/badge/code%20style-ruff-000000.svg)](https://github.com/astral-sh/ruff)

> **Production-grade ChIP-seq / ATAC-seq / Bisulfite-seq differential analysis framework
> built around the [ChIP-Atlas Differential Analysis API](https://chip-atlas.org).**

---

## Overview

`chipatlas-diff-pipeline` is an end-to-end Python framework for performing
differential chromatin binding (DiffBind) and differential DNA methylation (DMR)
analyses using public ChIP-seq, ATAC-seq, DNase-seq, or Bisulfite-seq data
curated in the [ChIP-Atlas](https://chip-atlas.org) database.

The pipeline:

1. **Submits** two groups of SRX/GSM experiment IDs to the ChIP-Atlas WABI API
2. **Polls** job status until completion (edgeR DPR or metilene DMR)
3. **Downloads** and extracts the ZIP result archive
4. **Parses** `wabi_result.bed` into an annotated pandas DataFrame
5. **Runs QC checks** — sex-chromosome confounds, CNV clusters, chrM contamination, MA asymmetry, low sample-size warnings
6. **Filters** non-standard contigs and implausibly small regions
7. **Generates** publication-quality volcano, MA, chromosome distribution, and size plots (plotnine + Prism theme)
8. **Annotates** top-50 regions with nearest gene symbols via the UCSC REST API
9. **Exports** CSVs, pickle object, and a Markdown summary report

### Why ChIP-Atlas?

ChIP-Atlas hosts uniformly processed ChIP-seq/ATAC-seq data from tens of thousands
of public experiments (SRA/GEO).  Their differential analysis API runs server-side
peak calling and statistical testing (edgeR for ChIP/ATAC, metilene for bisulfite),
making it possible to compare any two groups of public experiments without downloading
or processing raw sequencing data locally.

---

## Architecture

```
chipatlas-diff-pipeline/
├── src/
│   ├── api_client.py       # ChIP-Atlas WABI API: submit, poll, retrieve
│   ├── diff_analysis.py    # End-to-end workflow orchestrator
│   ├── parser.py           # wabi_result.bed → annotated pandas DataFrame
│   ├── filters.py          # QC & user-defined post-hoc filters
│   ├── qc_checks.py        # Automated QC suite (sex-chrom, CNV, chrM, MA…)
│   ├── visualization.py    # Volcano / MA / chromosome / size plots (plotnine)
│   ├── gene_annotation.py  # Nearest-gene lookup via UCSC REST API
│   ├── report.py           # Markdown report + CSV/pickle export
│   ├── cache.py            # Disk-backed analysis cache (SHA-256 keyed)
│   ├── validation.py       # Genome/ID/parameter validation
│   └── utils.py            # Logging, ID loading, directory helpers
├── configs/
│   ├── example_diffbind.yaml
│   └── example_dmr.yaml
├── examples/
│   ├── run_tp53_k562.py
│   ├── group_a_k562_dmso.txt
│   ├── group_b_k562_dauno.txt
│   └── sample_metadata.csv
├── tests/
│   └── test_pipeline.py
├── results/                  # Default output directory (git-ignored)
├── cli.py                    # Command-line interface
├── Dockerfile
├── docker-compose.yml
├── requirements.txt
├── pyproject.toml
└── README.md
```

---

## Installation

### Prerequisites

- Python 3.9 or newer
- Internet access (ChIP-Atlas API + optional UCSC gene annotation)

### From source

```bash
git clone https://github.com/your-org/chipatlas-diff-pipeline.git
cd chipatlas-diff-pipeline

# Create and activate virtual environment
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Install the CLI in editable mode
pip install -e .
```

### Verify installation

```bash
chipatlas-diff --help
python -c "from src.diff_analysis import run_diff_workflow; print('OK')"
```

---

## Quick Start

### Option 1 — Built-in example dataset

```bash
chipatlas-diff example --dataset tp53_k562_dmso_vs_daunorubicin
```

This submits a real TP53 ChIP-seq comparison (K562 DMSO vs Daunorubicin,
8 vs 8 experiments from GSE131484) and saves results to
`results/tp53_k562_dmso_vs_daunorubicin/`.

### Option 2 — YAML configuration file

```bash
chipatlas-diff run --config configs/example_diffbind.yaml
```

### Option 3 — Python API

```python
from src.diff_analysis import run_diff_workflow

results = run_diff_workflow(
    experiments_a = ["SRX5865959", "SRX5865960", "SRX5865961", "SRX5865962"],
    experiments_b = ["SRX5865967", "SRX5865968", "SRX5865969", "SRX5865970"],
    genome        = "hg38",
    analysis_type = "diffbind",
    description_a = "K562_DMSO",
    description_b = "K562_Daunorubicin",
    output_dir    = "results/my_analysis",
)

df  = results["diff_regions"]
sig = df[df["significant"]]
print(f"Significant differential regions: {len(sig)}")
```

### Option 4 — Two text files (CLI)

```bash
chipatlas-diff run \
    --group-a examples/group_a_k562_dmso.txt \
    --group-b examples/group_b_k562_dauno.txt \
    --genome hg38 \
    --analysis-type diffbind \
    --description-a K562_DMSO \
    --description-b K562_Daunorubicin \
    --output-dir results/my_analysis
```

---

## API Workflow

The ChIP-Atlas WABI API endpoint accepts a POST request with experiment IDs
and returns a request ID.  The pipeline then polls the status endpoint and
downloads the ZIP archive once the job finishes.

```
POST https://dtn1.ddbj.nig.ac.jp/wabi/chipatlas/
  ├── genome=hg38
  ├── antigenClass=diffbind
  ├── bedAFile=SRX001\nSRX002\n...
  └── bedBFile=SRX003\nSRX004\n...

→ requestId: wabi_20240101-123456-001-12345678-aa

GET .../wabi_20240101-...-aa?info=status
→ status: finished

GET .../wabi_20240101-...-aa?info=result&format=zip
→ ZIP containing wabi_result.bed, wabi_result.igv.bed, wabi_result.log
```

### BED output format

`wabi_result.bed` is a tab-separated file with 8 columns (no header):

| Column     | Description                                                              |
|------------|--------------------------------------------------------------------------|
| chrom      | Chromosome                                                               |
| chromStart | Region start (0-based)                                                   |
| chromEnd   | Region end                                                               |
| counts_A   | Comma-separated normalised read counts per Group A experiment            |
| counts_B   | Comma-separated normalised read counts per Group B experiment            |
| logFC      | Log₂ fold change (positive = A enriched, negative = B enriched)         |
| pvalue     | edgeR / metilene p-value                                                 |
| qvalue     | Benjamini-Hochberg FDR                                                   |

The pipeline adds: `region_size`, `direction`, `significant`, `mean_count_a`,
`mean_count_b`, and `score` (BED-compatible −log₁₀(q) × 100).

---

## Real Dataset Examples

### TP53 K562: DMSO vs Daunorubicin (primary — n=8 vs 8)

```yaml
# configs/example_diffbind.yaml
genome: hg38
analysis_type: diffbind
description_a: K562_DMSO
description_b: K562_Daunorubicin

experiments_a:  # K562 + DMSO
  - SRX5865959  - SRX5865960  - SRX5865961  - SRX5865962
  - SRX5865963  - SRX5865964  - SRX5865965  - SRX5865966

experiments_b:  # K562 + Daunorubicin
  - SRX5865967  - SRX5865968  - SRX5865969  - SRX5865970
  - SRX5865971  - SRX5865972  - SRX5865973  - SRX5865974
```

**Source:** GSE131484 (Bhatt et al.)  
**Biology:** Identifies genomic regions where daunorubicin-induced DNA damage
alters TP53 occupancy across 8 TP53 genotype backgrounds (WT, KO, 6 hotspot
mutants).  
**Expected runtime:** ~5–10 min  
**Design caveat:** Each group contains 8 different TP53 genotypes, not 8
replicates of the same genotype.

### TP53 MOLM-13 vs K-562 (cross-cell-line, QC demo — n=4 vs 4)

```python
experiments_a = ["SRX5865975", "SRX5865976", "SRX5865977", "SRX5865982"]  # MOLM-13 (male, AML)
experiments_b = ["SRX5865959", "SRX5865960", "SRX5865961", "SRX5865966"]  # K-562 (female, CML)
```

**Biology:** Cross-cell-line TP53 comparison under DMSO (matched genotypes).  
**QC demo:** Automatically detects and flags the expected chrY sex-chromosome
confound (MOLM-13 = male, K-562 = female).

---

## Supported Genomes

| Assembly  | Organism            | Standard chromosomes         |
|-----------|---------------------|------------------------------|
| `hg38`    | Human (GRCh38)      | chr1–22, X, Y, M             |
| `hg19`    | Human (GRCh37)      | chr1–22, X, Y, M             |
| `mm10`    | Mouse (GRCm38)      | chr1–19, X, Y, M             |
| `mm9`     | Mouse (NCBI37)      | chr1–19, X, Y, M             |
| `rn6`     | Rat (Rnor_6.0)      | chr1–20, X, Y, M             |
| `dm6`     | Drosophila (dm6)    | 2L, 2R, 3L, 3R, 4, X, Y, M  |
| `dm3`     | Drosophila (dm3)    | 2L, 2R, 3L, 3R, 4, X, Y, M  |
| `ce11`    | C. elegans (WBcel235)| I–V, X, M                   |
| `ce10`    | C. elegans (WS220)  | I–V, X, M                   |
| `sacCer3` | S. cerevisiae (R64) | I–XVI, M                     |

---

## Output Files

After a successful run, `output_dir/` contains:

```
output_dir/
├── analysis_object.pkl          # Complete results dict (downstream Python use)
├── diff_regions_all.csv         # All QC-filtered differential regions
├── diff_regions_significant.csv # FDR < 0.05 regions
├── diff_regions_top50.csv       # Top 50 by significance (with gene annotations)
├── diff_regions_unfiltered.csv  # Pre-QC regions (non-standard contigs included)
├── summary_report.md            # Human-readable Markdown report
├── volcano_plot.png / .svg      # logFC vs −log₁₀(Q-value)
├── chromosome_distribution.png / .svg
├── region_size_distribution.png / .svg
├── ma_plot.png / .svg           # log₂ mean count vs logFC
└── raw_results/
    ├── wabi_result.bed          # Original 8-column BED from ChIP-Atlas
    ├── wabi_result.igv.bed      # Colour-coded IGV visualisation BED
    ├── wabi_result.log          # edgeR / metilene analysis log
    └── wabi_result.igv.xml      # IGV session file
```

### Loading results in Python

```python
import pickle

with open("results/my_analysis/analysis_object.pkl", "rb") as fh:
    obj = pickle.load(fh)

df  = obj["diff_regions"]        # pd.DataFrame of all QC-filtered regions
sig = df[df["significant"]]      # FDR < 0.05 subset
print(obj["qc_warnings"])        # list of QC warning dicts
```

---

## Automated QC Checks

The pipeline runs eight automated QC checks on every result:

| Check | Severity | Description |
|-------|----------|-------------|
| Sex-chromosome confound (all regions) | **Major** | >80% one-directional chrX/Y with near-zero depleted counts → sex mismatch |
| Sex-chromosome confound (sig. regions) | **Major** | Same check restricted to significant regions |
| Non-standard contigs | Minor | Regions on unplaced scaffolds / random contigs |
| Implausibly small regions | Minor | < 10 bp — biologically implausible for TF binding |
| chrM enrichment | Major/Minor | >5 significant mitochondrial regions (ChIP/ATAC only) |
| Genomic clusters | Minor | ≥10 significant regions within 2 Mb → possible CNV |
| MA-plot asymmetry | Minor | Significant hits in one direction clustering at low counts |
| Region size distribution | Minor | Median significant region < 100 bp |
| Low sample size | Minor/Major | n < 5 or n ≤ 3 per group |

QC warnings are printed during the run and written in detail to `summary_report.md`.

---

## Visualisations

All four plots are exported as PNG (300 DPI) and SVG:

| Plot | Description |
|------|-------------|
| `volcano_plot` | logFC vs −log₁₀(Q-value); coloured by enrichment direction |
| `chromosome_distribution` | Grouped bar chart of significant regions per standard chromosome |
| `region_size_distribution` | Log-scale histogram of region sizes by direction |
| `ma_plot` | log₂ mean normalised count (A) vs logFC — reveals low-count noise |

Colours follow a consistent palette: **red** = A-enriched, **blue** = B-enriched,
**grey** = not significant.

---

## CLI Reference

```
chipatlas-diff <command> [options]

Commands:
  run       Submit and run a differential analysis
  example   Run a built-in example dataset
  cache     Manage the local analysis cache

chipatlas-diff run [options]:
  --config FILE             YAML configuration file
  --group-a FILE            Group A IDs (one per line)
  --group-b FILE            Group B IDs (one per line)
  --group-file FILE         Two-group metadata CSV/TSV
  --genome GENOME           Assembly (default: hg38)
  --analysis-type TYPE      diffbind or dmr (default: diffbind)
  --title TITLE             Short analysis title
  --description-a LABEL     Group A label (default: group_A)
  --description-b LABEL     Group B label (default: group_B)
  --output-dir DIR          Output directory (default: diff_analysis_results)
  --min-score INT           Minimum BED score 0–1000 (default: 0)
  --min-region-size BP      Minimum region size in bp (default: 0)
  --direction-filter DIR    A_enriched or B_enriched (default: both)
  --no-cache                Force API re-submission (ignore cache)
  --no-plots                Skip visualisation generation
  --no-genes                Skip UCSC gene annotation
  --job-timeout SEC         Max wait time in seconds (default: 900)
  --poll-interval SEC       Status poll interval (default: 15)
  --log-level LEVEL         DEBUG/INFO/WARNING/ERROR (default: INFO)

chipatlas-diff example [options]:
  --dataset NAME            tp53_k562_dmso_vs_daunorubicin (default)
                            tp53_molm13_vs_k562
  --output-dir DIR          Override default output path

chipatlas-diff cache [list|clear]:
  list                      List all cached analysis keys
  clear                     Delete all cached analyses
```

---

## Docker Usage

### Build the image

```bash
docker build -t chipatlas-diff-pipeline:latest .
```

### Run a built-in example

```bash
docker run --rm \
    -v "$(pwd)/results:/app/results" \
    chipatlas-diff-pipeline:latest \
    example --dataset tp53_k562_dmso_vs_daunorubicin
```

### Run from a YAML config

```bash
docker run --rm \
    -v "$(pwd)/results:/app/results" \
    -v "$(pwd)/configs:/app/configs:ro" \
    chipatlas-diff-pipeline:latest \
    run --config /app/configs/example_diffbind.yaml
```

### Using docker compose

```bash
# Start analysis service
docker compose run analysis --config /app/configs/example_diffbind.yaml

# List cached analyses
docker compose run chipatlas-diff cache list
```

---

## Running Tests

```bash
# Install dev dependencies
pip install -e ".[dev]"

# Run all tests
pytest tests/ -v

# Run with coverage report
pytest tests/ -v --cov=src --cov-report=term-missing

# Run a specific test module
pytest tests/test_pipeline.py::TestParser -v
```

The test suite covers:
- Input validation (genome names, ID formats, minimum replicates)
- BED parsing (8-column, IGV BED9, edge cases)
- Filtering (standard chromosomes, min size, FDR, logFC, direction)
- QC checks (sex-chromosome detection, cluster detection, sample size)
- Cache (set/get/delete/list, deterministic keys)
- Utility functions (file loading, two-group CSV)

---

## Performance Notes

| Scenario | Approximate runtime |
|----------|---------------------|
| 2 vs 2 experiments | 2–3 min |
| 4 vs 4 experiments | 3–5 min |
| 8 vs 8 experiments | 5–10 min |
| 16 vs 16 experiments | 10–20 min |

- Job runtimes are dominated by server-side peak calling and edgeR fitting,
  not network transfer.
- The local cache (SHA-256 keyed by experiment IDs + genome + analysis type)
  means identical analyses run instantly on repeated calls.
- UCSC gene annotation queries the REST API for up to 50 regions sequentially
  (≈ 1–5 s total). Disable with `--no-genes` if speed is critical.

---

## Caching

Completed analyses are cached under `~/.chipatlas_cache/` using a 16-character
hex key derived from the sorted experiment IDs, genome, and analysis type.
This prevents duplicate API submissions for identical inputs.

```python
# Clear cache programmatically
from src.cache import AnalysisCache
cache = AnalysisCache()
cache.clear()

# Or via CLI
chipatlas-diff cache clear
```

Pass `use_cache=False` (Python) or `--no-cache` (CLI) to force a new submission.

---

## Troubleshooting

### API submission fails

```
RuntimeError: ChIP-Atlas API unreachable
```

- Check internet connectivity: `curl https://chip-atlas.org`
- Check API status at https://chip-atlas.org
- Some institutional networks block outbound HTTPS to DDBJ servers — try from
  a different network or use the Docker image with a VPN

### Job times out

```
RuntimeError: Job ... timed out after 900s
```

- Increase the timeout: `--job-timeout 1800`
- Large comparisons (≥16 vs 16) may need up to 30 min

### No BED file in ZIP

```
RuntimeError: No BED file found in ZIP archive
```

- The job may have failed server-side — check the `raw_results/wabi_result.log`
  if partial extraction occurred
- Experiment IDs not present in ChIP-Atlas (not all SRX IDs are processed)

### Experiment IDs not recognised

- Only SRX/ERX/DRX/GSM IDs from experiments processed by ChIP-Atlas are valid
- Search at https://chip-atlas.org to verify ID availability and coverage
- GSM IDs are automatically mapped to SRX by ChIP-Atlas

### All results on non-standard contigs

- Verify the `--genome` flag matches the assembly used in the original experiments
- Human experiments always require `hg38` or `hg19`

### Sex-chromosome warning on human data

- Expected when comparing male vs female cell lines (e.g. HeLa [female] vs
  HEK293 [female, but chrY present]; K562 [female] vs MOLM-13 [male])
- Filter chrX/chrY from the output or note the confound in your analysis

---

## Python API Reference

### `run_diff_workflow`

```python
from src.diff_analysis import run_diff_workflow

results = run_diff_workflow(
    experiments_a   = [...],      # list of str, ≥ 2 SRX/GSM IDs
    experiments_b   = [...],      # list of str, ≥ 2 SRX/GSM IDs
    genome          = "hg38",     # see Supported Genomes table
    analysis_type   = "diffbind", # or "dmr"
    title           = "my_run",
    description_a   = "Condition_A",
    description_b   = "Condition_B",
    min_score       = 0,          # post-filter: BED score ≥ N
    min_region_size = 0,          # post-filter: size ≥ N bp
    direction_filter = None,      # "A_enriched" | "B_enriched" | None
    design_caveats  = [],         # list of str, appended to report
    output_dir      = "results/", # output directory
    use_cache       = True,       # check cache before submitting
    generate_plots  = True,       # generate plotnine visualisations
    annotate_genes  = True,       # query UCSC for gene symbols
    job_timeout     = 900,        # seconds
    poll_interval   = 15,         # seconds
)
```

Returns a dict with keys:
- `diff_regions` — `pd.DataFrame` (QC-filtered)
- `diff_regions_unfiltered` — `pd.DataFrame` (pre-QC)
- `qc_warnings` — `list[dict]`
- `experiments_a`, `experiments_b` — input accession lists
- `raw_files` — dict of extracted file paths
- `log_content` — edgeR/metilene log string
- `parameters` — echo of all parameters

---

## Citation

If you use this pipeline in a publication, please cite:

**ChIP-Atlas:**
> Zou Z, Ohta T, Oki S. (2024). ChIP-Atlas 3.0: a data-mining suite to
> explore chromatin accessibility and chromosome architecture together with
> histone modifications. *Nucleic Acids Research*, 52(W1):W45–W53.
> https://doi.org/10.1093/nar/gkae358

**edgeR (for diffbind results):**
> Robinson MD, McCarthy DJ, Smyth GK. (2010). edgeR: a Bioconductor package
> for differential expression analysis of digital gene expression data.
> *Bioinformatics*, 26(1):139–140.
> https://doi.org/10.1093/bioinformatics/btp616

**metilene (for DMR results):**
> Jühling F, et al. (2016). metilene: Fast and sensitive calling of
> differentially methylated regions from bisulfite sequencing data.
> *Genome Research*, 26(2):256–262.
> https://doi.org/10.1101/gr.196394.115

---

## Contributing

Contributions are welcome!

```bash
# Fork and clone
git clone https://github.com/your-org/chipatlas-diff-pipeline.git
cd chipatlas-diff-pipeline

# Create a feature branch
git checkout -b feat/my-feature

# Install dev dependencies
pip install -e ".[dev]"

# Run linter and tests
ruff check src/ cli.py
pytest tests/ -v

# Submit a pull request
```

Please follow the existing code style (typed functions, NumPy-format docstrings,
descriptive variable names).

---

## License

This project is released under the [MIT License](LICENSE).

ChIP-Atlas data is made available under the terms described at
https://chip-atlas.org/data_download.  Users are responsible for ensuring
their use of public ChIP-seq data complies with the original data policies.
