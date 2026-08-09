# mtDNA-ML-Predictor

A machine-learning pipeline for predicting low-frequency mitochondrial DNA (mtDNA) mutations from [mtDNApipe](https://github.com/FMMU-Xing-Lab/mtDNApipe) variant-calling outputs.

The pipeline takes the per-sample mitochondrial BAM and mutation files (hetro/homo) produced by mtDNApipe, extracts 29 alignment/sequencing features, removes 6 features that are no longer used, adds 9 annotation features, and finally trains an XGBoost model on the HCC tissue training set to predict every candidate mutation in new samples.

## Workflow

```mermaid
flowchart LR
    A[mtDNApipe output<br/>per sample: mt.no.softclip.bam + bai<br/>+ hetro/homo mutation files] --> B[1. Feature extraction<br/>get_output_true_file.R<br/>+ get_feature_mtDNApipe.py<br/>+ get_total_output_features.R]
    B --> C[2. Remove 6 features<br/>AF, dp, mosaic/het/refhom/althom likelihood]
    C --> D[3. Add 9 annotation features<br/>repeat-region, Population-freq, NAV, type,<br/>if_Trans, if_version, region, VAF_mitomap, dbSNP]
    D --> E[4. XGBoost prediction<br/>training set: HCC_training_all.txt]
    E --> F[Output: mutation + pred_label + pred_prob]
```

The final prediction output contains only three columns: the mutation ID (`sample~chrM~position~reference_base~mutant_base`), the predicted class (0/1), and the probability of being a real mutation.

## Directory Layout

```text
mtDNA-ML-Predictor/
├── run_pipeline.sh            # One-command master pipeline (recommended entry point)
├── check_environment.sh       # Environment check (runs automatically before the pipeline)
├── download_reference.sh      # Re-download the mtDNA reference genome (rCRS)
├── requirements.txt           # pip dependencies
├── environment.yml            # conda environment (recommended)
├── LICENSE
├── README.md
├── scripts/
│   ├── organize_samples.sh          # Organize sample directories / create sample_name.txt
│   ├── get_feature_pipeline.sh      # Feature extraction (portable version of get_mtDNApipe_feature.sh)
│   ├── get_output_true_file.R       # Prepare candidate variant sites output_true_<vaf>
│   ├── get_feature_mtDNApipe.py     # Extract per-sample alignment features (parallel)
│   ├── get_total_output_features.R  # Aggregate features across all samples
│   ├── delete_features.py           # Remove the 6 unused features (replaces the old manual step)
│   ├── add_features.py              # Add 9 annotation features (portable version of 添加10个特征.py)
│   └── predict.py                   # Train + predict (portable version of 预测集验证.py, trimmed output)
├── data/
│   ├── HCC_training_all.txt         # HCC training set (formerly 肝癌训练集-all.txt)
│   ├── dbSNP.txt
│   ├── mitomap.txt
│   ├── mitomap-snp.txt
│   ├── region.txt
│   └── mtDNA_region_ge5.txt         # formerly mtDNA区间≥5.txt
└── reference/
    └── human_mtDNA.fasta            # rCRS (NC_012920.1), sequence name chrM, 16,569 bp
```

## Requirements & Installation

A Linux server is required. The following tools must be available:

- `python3` (>= 3.7) with: `pysam`, `pandas`, `numpy`, `scipy`, `regex`, `pyfaidx`, `xgboost`, `scikit-learn`
- `Rscript` (the R scripts only use base R; no additional R packages are needed)
- `samtools` (builds the reference index and computes read lengths)

Recommended: create the conda environment in one step:

```bash
conda env create -f environment.yml
conda activate mtdna-ml
```

Or install with pip:

```bash
python3 -m pip install -r requirements.txt
# samtools and R can be installed with your system package manager, e.g.:
# sudo apt-get install -y samtools r-base
```

Before running, you can automatically check whether the current environment is ready:

```bash
bash check_environment.sh
```

The master pipeline `run_pipeline.sh` runs this environment check automatically by default (skip it with `-k`).

## Input Requirements

Put all samples in one analysis folder (e.g., `sample_dir/`). Either of the following layouts is accepted:

Layout 1 (recommended: one sub-folder per sample):

```text
sample_dir/
├── S1/
│   ├── S1.mt.no.softclip.bam
│   ├── S1.mt.no.softclip.bam.bai   (or S1.bai)
│   ├── S1.hetro_0.1.txt
│   └── S1.homo_0.1.txt
└── S2/
    └── ...
```

Layout 2 (flat; the pipeline will organize it into Layout 1 automatically):

```text
sample_dir/
├── S1.mt.no.softclip.bam
├── S1.bai
├── S1.hetro_0.1.txt
├── S1.homo_0.1.txt
└── ...
```

Notes:

- Samples are identified by `*.mt.no.softclip.bam` files; `*.mt.bam` names are also accepted automatically.
- The `0.1` in the hetro/homo file names is the VAF tag and must match the `-v` argument (default: `0.1`).
- The reference genome must be the same reference used for alignment (the bundled default is rCRS with the sequence named `chrM`).

## Quick Start

```bash
bash run_pipeline.sh -i /path/to/sample_dir -t 8
```

Common options:

```bash
bash run_pipeline.sh \
    -i /path/to/sample_dir \          # analysis folder containing the samples
    -o /path/to/outputs \             # output directory for final results (default: sample_dir/outputs)
    -c 2 \                            # cutoff (%) for low-frequency variants; no filtering by default,
                                      # enable it with the environment variable APPLY_CUTOFF=1
    -v 0.1 \                          # VAF tag used in the hetro/homo file names
    -t 8 \                            # number of parallel feature-extraction jobs
    -r /path/to/human_mtDNA.fasta \   # reference genome (default: bundled rCRS)
    -T /path/to/HCC_training_all.txt \# training set (default: bundled data)
    -s data1                        # sample/cohort name (used in the log)
```

Run the steps individually (equivalent to the master pipeline):

```bash
# 1) Feature extraction
bash scripts/get_feature_pipeline.sh /path/to/sample_dir 2 8 0.1 reference/human_mtDNA.fasta

# 2) Remove the 6 unused features
python3 scripts/delete_features.py \
    --input /path/to/sample_dir/Total_output_feature_0.1.txt \
    --output outputs/feature.txt

# 3) Add the annotation features
python3 scripts/add_features.py \
    --input outputs/feature.txt \
    --output outputs/feature_add10.txt \
    --data-dir data

# 4) Train + predict
python3 scripts/predict.py \
    --train data/HCC_training_all.txt \
    --predict outputs/feature_add10.txt \
    --output outputs/data1-pred-results.txt \
    --sample-name data1
```

## Output Files

| File | Description |
| --- | --- |
| `<sample_dir>/output_true_<vaf>` | Candidate variant sites per sample (chrM, position, bases, frequency) |
| `<sample_dir>/Total_output_feature_<vaf>.txt` | Aggregated 29-column raw features for all samples |
| `<output>/feature.txt` | 23-column feature table after removing the 6 features |
| `<output>/feature_add10.txt` | 32-column feature table after adding the 9 annotation features |
| `<output>/<sample_set>-pred-results.txt` | Final predictions: `id`, `pred_label`, `pred_prob` |

### Feature Columns

Of the 29 features extracted per variant, the following 6 are removed before training (the default columns dropped by `delete_features.py`; customize with `--drop`):

```text
AF, dp, mosaic_likelihood, het_likelihood, refhom_likelihood, althom_likelihood
```

The 9 annotation features added by `add_features.py`:

```text
repeat-region, Population-freq, NAV, type, if_Trans, if_version, region, VAF_mitomap, dbSNP
```

The final 32 columns exactly match `HCC_training_all.txt` (after removing the `label` column).

## FAQ

**1. The environment check fails. What should I do?**

Install the missing command or Python packages as indicated: `pip install -r requirements.txt`; install `samtools` and `R` with `apt-get`/`conda`.

**2. The reference index is missing?**

The master pipeline automatically creates it with `samtools faidx reference/human_mtDNA.fasta` (generates the `.fai` file).

**3. My BAM is not named `*.mt.no.softclip.bam`?**

The scripts also try to match `*.mt.bam`. If no BAM is found, check `sample_name.txt` and the BAM file names.

**4. How can I restore or change the deleted features?**

`delete_features.py` accepts a custom `--drop` list, and the annotation features added by `add_features.py` can also be modified as needed.

**5. Mapping to the original scripts**

| Original file | This repository |
| --- | --- |
| `get_mtDNApipe_feature.sh` | `scripts/get_feature_pipeline.sh` |
| `get_output_true_file.R` | `scripts/get_output_true_file.R` |
| `get_feature_mtDNApipe.py` | `scripts/get_feature_mtDNApipe.py` (reference path parameterized) |
| `get_total_output_features.R` | `scripts/get_total_output_features.R` |
| `添加10个特征.py` | `scripts/add_features.py` |
| `预测集验证.py` | `scripts/predict.py` (only mutation IDs and predictions are output) |
| `肝癌训练集-all.txt` | `data/HCC_training_all.txt` |
| `mtDNA区间≥5.txt` | `data/mtDNA_region_ge5.txt` |

**6. Prediction results differ slightly from the old `Figure7-pred-results.txt`?**

The feature-alignment logic is identical, but different xgboost versions introduce small floating-point differences during training. The differences are concentrated on mutations whose predicted probability is close to 0.5 (in this repository's validation: only ~1.3% of the 304k mutations flipped at the 0.5 boundary, with a mean probability difference of ~0.018). To reproduce the old results exactly, install the same xgboost version used in the original analysis.

## Citation

If you use this pipeline, please cite mtDNApipe and the original paper of this project (to be added once published).
