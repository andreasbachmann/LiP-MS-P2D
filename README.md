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
git clone https://github.com/gruber-sciencelab/LiP-MS-P2D.git
cd LiP-MS-P2D

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
Rscript scripts/msstatslip.R --test

# Run domain-level aggregation on test data
Rscript scripts/combine.R --test
```

Test outputs will be in `test/results/`.

# Full Analysis

## Required Input Files

1. **Spectronaut output**: LiP and TrP fragment-level quantification
2. **FASTA file**: Protein sequences for your organism (same one as used in Spectronaut)
3. **Domain annotations**: From InterPro (see below)

## Step 1: MSstatsLiP Processing

```bash
# Edit input paths in scripts/msstatslip.R, then run:
Rscript scripts/msstatslip.R
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
- `results/domain_level_results_ebm_dedup.csv`: Deduplicated EBM results
- `results/domain_level_results_cauchy_dedup.csv`: Deduplicated Cauchy results
- `results/domain_level_results_ebm_dedup_global.csv`: Global deduplicated EBM results
- `results/domain_level_results_cauchy_dedup_global.csv`: Global deduplicated Cauchy results
- `results/peptide_domain_mapping.csv`: Peptide-to-domain mapping details

# Custom Domain Annotation

A predefined domain annotation file for human proteins is provided for the analysis (in `data/domain_annotations_consolidated.csv`). However, if you work with another model organism, you can create your own custom domain annotation file as described below.

## Step 1: Download Required Files

**InterPro protein annotations:**
Download the `protein2ipr.dat.gz` from the InterPro FTP server and save it under `data/`.
```bash
wget https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/protein2ipr.dat.gz
gunzip protein2ipr.dat.gz
```
Beware that this file contains domain annotations for all UniProtKB proteins and will be about 100GB unzipped. 

**Organism-specific protein list:**
Download your organism's proteome from UniProtKB. For example, for human (UP000005640), navigate to `https://www.uniprot.org/proteomes/UP000005640`, under Components, click on Download. Select .tsv as format and deselect all columns (under Customize columns) except for Entry Name. Then download and save that file under `data/`.

## Step 2: Subset to Your Organism

Filter the InterPro annotations to keep only proteins from your organism:
```bash
awk -F '\t' 'NR==FNR { proteins[$1]; next } $1 in proteins' \
    data/your_organism_proteins.list data/protein2ipr.dat > data/protein_annotation_organism.tsv
```

## Step 3: Consolidate Overlapping Domains

The raw InterPro file contains redundant entries because the same domain region is often annotated by multiple source databases (Pfam, SMART, PROSITE, etc.) with slightly different names/boundaries. The `cluster_domains.R` script merges these overlapping annotations. Your first argument specifies the input, while the second argument specifies the output.

```bash
Rscript scripts/cluster_domains.R \
    data/protein_annotation_organism.tsv \
    data/domain_annotations_consolidated.csv
```

The output file `domain_annotations_consolidated.csv` can be used directly as input for `combine.R`.


# RNA-binding domain Compilation

To test if domains are enriched in RNA-binding, a list of RNA-binding domains has been compiled from multiple sources.

## Input Files

| File | Description | Source |
|------|-------------|--------|
| `protein_annotation_human.tsv` | Human protein domain annotations from InterPro | Generated in the domain annotation step |
| `interpro2go` | InterPro to Gene Ontology mapping | [InterPro FTP](https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/) |
| `Search_Results_InterPro.csv` | RNA-binding domain annotations | Exported from [RBP2GO](https://rbp2gov2.dkfz.de/) |

## Sources

1. **Keyword search**: The human domain annotations were searched for entries matching RNA-related terms (e.g. 'RNA recognition motif', 'RNA-binding', 'RRM', 'DNA/RNA-binding'). Domains belonging to the InterPro RNA binding superfamily (IPR035979) were also manually included.

2. **InterPro2GO mapping**: The InterPro2GO file was downloaded from the [InterPro FTP server](https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/). All InterPro entries annotated with the GO term 'RNA binding' were extracted.

3. **RBP2GO database**: Curated RNA-binding domain identifiers were obtained from [RBP2GO](https://rbp2go.dkfz.de/) (Herger et al., 2021) and combined with the InterPro-based annotations.

## Usage
```bash
# Step 1: Compile InterPro-based RNA-binding domain IDs
bash scripts/compile_RNAbindingIDs.sh

# Step 2: Combine with RBP2GO annotations and deduplicate
Rscript scripts/combine_RNAs.R
```

The final deduplicated list (`unique_rbd.tsv`) is used to check if a domain hit is an RNA-binding domain for the Fisher's exact test.