# Rat TWAS Hub

Code to generate the data and website for Rat TWAS Hub

Code for running large-scale multi-modal TWAS analyses in rat and producing static Markdown reports. This repository contains scripts that compute TWAS statistics for many traits, merge and annotate them, produce reports, and generate an interactive website. This code is adapted from the (human) [TWAS HUB](https://github.com/gusevlab/TWAS_HUB). It extends gene expression-based TWAS to multimodal RNA phenotypes (xTWAS) using [Pantry](https://github.com/PejLab/Pantry) to generate the transcriptome models.

## Parameter files

The code uses the following parameter files:
 
* `panels.par`: Statistics on the transcriptome reference panels used. A `MODALITY` column has been added for xTWAS.
* `all.models.par`: Statistics on each of the predictive models. The `Snakefile` assembles this from the FUSION profiler script outputs.
* `traits.par`: Information on each of the traits/TWAS studies performed. Important, the `OUTPUT` column must to point to the merged results from all chromosomes for that trait.
* `gene_names.tsv`: Lookup table to convert the Ensembl IDs in the results to gene symbols for the site. The `Snakefile` assembles this from a gene annotation GTF file.

## Workflow

* Use Snakemake with the included `Snakefile` to run TWAS ([FUSION](https://github.com/gusevlab/fusion_twas)) and post-processing.
* Run `scripts/built_jekyll.sh` to generate the Jekyll template.
* `cd jekyll` and then `bundle exec jekyll build` to generate the site.
