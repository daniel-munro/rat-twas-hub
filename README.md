# Rat TWAS Hub

Code to generate the data and website for Rat TWAS Hub.

This repository runs large-scale, multi-modal TWAS analyses in rat and produces static Markdown reports plus a Jekyll site. It is adapted from the (human) [TWAS HUB](https://github.com/gusevlab/TWAS_HUB) and extends gene expression-based TWAS to multi-modal RNA phenotypes (xTWAS) using [Pantry](https://github.com/PejLab/Pantry) to generate phenotypes for the transcriptome models.

## Requirements

- Python 3 with `snakemake` and `pandas`.
- R (4.x recommended) with packages:
  - required: `tidyverse`, `yaml`, `biomaRt`, `optparse`, `plink2R`, `RColorBrewer`, `coloc`
- GNU core tools (`awk`, `sed`, `join`, `tar`, `gzip`).
- Ruby + `bundler` for Jekyll site builds (see `jekyll/Gemfile`).
- `rsync` + `ssh` for deployment (optional).

Notes:

- Fetching the cross-species ortholog mapping uses Ensembl via `biomaRt` and requires network access. Routine site updates use the saved mapping.
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

- `data/all_models.tsv`, `data/twas_out/`, `data/traits.summary.tsv`, `data/genes_n_models.tsv`, `data/genes_n_assoc.tsv`, `data/site_loci.tsv`
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

To generate only TWAS result tables, without site outputs, use `twas_results`:

```bash
snakemake --profile slurm --dry-run --printshellcmds twas_results
snakemake --profile slurm twas_results
```

This produces `data/twas_out/{trait}.{all,top,post}.tsv` for all traits in
`data/traits.tsv`. Add `--config project=PROJECT_ID` to select one project.
The default target builds results and site source data for all traits;
project selection applies only to `twas_results`.

Plots and download archives are generated per trait. To bundle short jobs into
fewer Slurm submissions, use [job grouping](https://snakemake.readthedocs.io/en/stable/executing/grouping.html)
with `--group-components`:

```bash
snakemake --profile slurm --dry-run --group-components site_archives=16 site_plots=8
snakemake --profile slurm --group-components site_archives=16 site_plots=8
```

Each submission contains up to 16 pending archive jobs or 8 pending plot jobs.
Only traits needing a rebuild are included. Group resource requests account for
the member jobs; each plot requests 2 GB and each archive 512 MB. Without
`--group-components`, each trait runs separately. Local execution ignores groups.

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

## Cross-species mapping

`data/cross_species/human_rat_orthologs.ensembl108.tsv` stores the ortholog
mapping separately from association counts. The initial frozen table was
extracted from the existing `jekyll/cross-species.json`, preserving its mapping
from the Ensembl 108 (October 2022) lookup. Keep this TSV with the analysis inputs.
Site updates join the current association counts onto this table without
contacting Ensembl.

If the mapping is absent, `cross_species_orthologs` fetches it. Model-list
modification times do not trigger a refresh; when intentionally changing the
model gene set or refreshing the mapping, save the old table and explicitly run:

```bash
snakemake --profile slurm --dry-run --forcerun cross_species_orthologs -- cross_species_orthologs
snakemake --profile slurm --forcerun cross_species_orthologs -- cross_species_orthologs
```

Trait-page hit counts and transcriptome-wide pleiotropy tests use
`0.05 / total model count` from `data/all_models.tsv`, including models with
missing TWAS results, matching TWAS hit selection and porcupine plots.

## Locus pages

`site_locus_manifest` is a local checkpoint: once all configured traits'
post-processing reports are available, it writes `data/site_loci.tsv` from their
rows. This manifest enumerates each trait's expected locus pages, conditional
CSVs, and model tables. Loci are numbered by chromosome and
start position, matching the site tables. No file globbing is used to discover
loci; a header-only report requires no locus pages.

Trait pages and their locus pages share a `site_trait_{trait}` job group, so a
trait rebuild bundles its locus rendering with the trait page. Locus files are
explicit outputs: a deleted page, CSV, or model table is rebuilt on the next
full workflow run. The trait-page rule removes that trait's generated locus
page and model-table directories before rebuilding, so obsolete loci are removed
when the report changes, including when the trait no longer has any loci.

The `.genes` and `.cond.csv` files named in the report are explicit inputs.
Missing source files fail the build; regenerate the affected chromosome with
`twas_post` if these original post-processing files have been removed.
