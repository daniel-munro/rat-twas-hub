# Rat TWAS Hub

Code to generate the data and website for Rat TWAS Hub.

This repository runs large-scale, multi-modal TWAS analyses in rat and produces static Markdown reports plus a Jekyll site. It is adapted from the (human) [TWAS HUB](https://github.com/gusevlab/TWAS_HUB) and extends gene expression-based TWAS to multi-modal RNA phenotypes (xTWAS) using [Pantry](https://github.com/PejLab/Pantry) to generate phenotypes for the transcriptome models.

## Requirements

- Python 3 with `snakemake` and `pandas`.
- R (4.x recommended) with packages:
  - required: `tidyverse`, `yaml`, `biomaRt`, `optparse`, `plink2R`, `RColorBrewer`, `coloc`
- GNU core tools (`awk`, `sed`, `join`, `tar`, `gzip`) and `parallel`.
- Ruby + `bundler` for Jekyll site builds (see `jekyll/Gemfile`).
- `rsync` + `ssh` for deployment (optional).

Notes:

- The cross-species step uses Ensembl via `biomaRt` and requires network access.
- TWAS runs require binary PLINK reference panels (see inputs below).

## Inputs and parameter files

The pipeline expects the following inputs under `data/`:

- `panels.tsv`: Transcriptome panel metadata (columns include `PANEL`, `TISSUE`, `MODALITY`, `N`).
- `traits.tsv`: Trait metadata (columns include `OUTPUT`, `ID`, `N`, `PROJECT`, `TAGS`, `NAME`, `DESCRIPTION`).
- `projects.tsv`: Project metadata used for the site.
- `sequence_report.tsv`: Chromosome sizes for porcupine plots.
- `WEIGHTS/<tissue>/<modality>.pos` and `<modality>.profile`: FUSION profiler outputs used to build `data/all_models.tsv`.
- `LDREF/Brain_v4.chr{chrom}.{bed,bim,fam}`: PLINK reference LD panels used by FUSION.
- `gwas_original/sumstats/{project}/regressedlr_{trait}_chrgwas{chrom}.mlma.gz`: GWAS summary statistics (GCTA output).
- `gwas_original/traits.tsv`: Project + trait metadata used by `scripts/prepare/prepare_trait_table.R` (columns: `project_id`, `trait_id`, `tags`, `name`, `description`).
- `cross_species/genes.models.nfo` and `cross_species/genes.nfo`: Human reference lists for the cross-species table.

Generated outputs live in:

- `data/all_models.tsv`, `data/twas_out/`, `data/traits.summary.tsv`, `data/genes_n_models.tsv`, `data/genes_n_assoc.tsv`
- `jekyll/` (traits, genes, and `_data` assets consumed by the site)

## Workflow

1. Prepare trait metadata (optional, if starting from raw GWAS):

```bash
scripts/prepare/prepare_traits.sh
```

This reads `data/gwas_original/traits.tsv` and GCTA logs to produce `data/traits.tsv`.

2. Run TWAS + post-processing + site data build:

```bash
snakemake -j 8
```

If running on HPC, be sure to run in a job, or have snakemake run steps as jobs using, e.g., [snakemake-executor-plugin-slurm](https://github.com/snakemake/snakemake-executor-plugin-slurm).

3. Build the Jekyll site:

```bash
cd jekyll
bundle install
bundle exec jekyll build
```

4. Deploy (optional, requires SSH access to the target host):

```bash
scripts/site/deploy.sh
```

## Repository layout

- `Snakefile`: Main workflow (TWAS, post-processing, and site data generation).
- `scripts/twas/`: FUSION TWAS and post-processing scripts.
- `scripts/site/`: Jekyll data builders and plot generation.
- `scripts/prepare/`: Trait table helpers for GWAS inputs.
- `jekyll/`: Site templates and assets.
