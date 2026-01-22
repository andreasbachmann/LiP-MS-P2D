# Overview

LiP-MS measures protein structural changes by detecting differences in protease accessibility. This pipeline:

1. Runs MSstatsLiP with Spectronaut data to generate peptide-level output
2. Filters to semi-tryptic peptides
3. Maps peptides to protein domains (using InterPro domain annotations)
4. Aggregates p-values using Cauchy combination test and Empirical Brown's Method (EBM)
5. Outputs domain-level structural changes

# Installation

## Prerequisites

- Conda (Miniconda or Anaconda)
- Git

## Setup

```bash
# Clone the repository and change into the project directory
git clone https://github.com/andreasbachmann/LiP-MS-P2D.git
cd lip-ms-p2d

# Create a new conda environment and activate it
conda create -n lip-ms-p2d -c conda-forge --file requirements.txt
conda activate lip-ms-p2d

# Install R packages not available via conda-forge
Rscript packages.R
```


# Validation Run

Verify the installation with test data:

```bash
# Run MSstatsLiP processing on test data
Rscript scripts/msstatslip_run.R --test

# Run domain-level aggregation on test data
Rscript scripts/combine.R --test
```

Test outputs will be in `test/results/`.

# Full Analysis

## Required Input Files

1. **Spectronaut output**: LiP and TrP fragment-level quantification
2. **FASTA file**: Protein sequences for your organism (same one as used in Spectronaut)
3. **Domain annotations**: From InterPro (see below)

## Step 1: MSstatsLiP Preprocessing

```bash
# Edit paths in scripts/msstatslip_run.R, then run:
Rscript scripts/msstatslip_run.R
```

**Inputs:**
- `data/lip_raw.csv`: Spectronaut LiP output
- `data/trp_raw.csv`: Spectronaut TrP output
- `data/protein.fasta`: Protein sequences

**Outputs:**
- `results/summarized_lip_data.csv`: Replicate intensity matrix (used later in EBM)
- `results/adjusted_lip_model.csv`: Peptide-level fold changes and p-values


## Step 2: Domain-Level Aggregation

```bash
Rscript scripts/combine.R
```

**Inputs:**
- `results/summarized_lip_data.csv`: From Step 1
- `results/adjusted_lip_model.csv`: From Step 1
- `data/domain_annotations_consolidated.csv`: InterPro domains

**Outputs:**
- `results/domain_level_results_ebm.csv`: EBM aggregation results
- `results/domain_level_results_cauchy.csv`: Cauchy aggregation results
- `results/domain_level_results_ebm_dedup.csv`: Deduplicated EBM results
- `results/domain_level_results_cauchy_dedup.csv`: Deduplicated Cauchy results
- `results/peptide_domain_mapping.csv`: Peptide-to-domain mapping details