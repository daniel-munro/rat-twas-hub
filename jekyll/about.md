---
title: About / FAQ
permalink: about/
---

# About

TWAS measures the genetic association between a transcriptomic phenotype such as gene expression and a complex phenotype using only GWAS summary-level data (see: [Gusev et al. 2016 *Nature Genetics*](https://www.ncbi.nlm.nih.gov/pubmed/26854917)). The TWAS central dogma is that associated genes are more likely to be causal mediators of the trait and thus informative of the trait's biological nature or as targets for experimental follow-up.

Rat TWAS Hub is an interactive browser of results from integrative analyses of GWAS and functional data.
The aim is facilitate the investigation of individual TWAS associations; pleiotropic trait associations for a given gene of interest; predicted gene associations for a given trait of interest with detailed per-locus statistics; and pleiotropic relationships between traits based on shared associated genes. See the [USAGE](/usage){: .border} tab for detailed examples of each analysis type.

Rat TWAS Hub is managed in the labs of [Abraham Palmer](https://palmerlab.org/) at UC San Diego and [Pejman Mohammadi](https://pejlab.org/) at Seattle Children's Research Institute and University of Washington. It was adapted from the (human) TWAS hub developed in the [Gusev Lab](http://gusevlab.org) at the Dana-Farber Cancer Institute and Harvard Medical School. For questions or comments about Rat TWAS Hub, please contact the developer and maintainer Daniel Munro at [dmunro@health.ucsd.edu](mailto:dmunro@health.ucsd.edu).

# FAQ

#### How is Rat TWAS hub generated?

For each trait, a TWAS is carried out using the [FUSION](http://gusevlab.org/projects/fusion/) software. FUSON post-processing is then used to extract all significant associations (after Bonferroni correction) and grouped into contiguous loci and a step-wise conditional analysis is performed to identify independent associations (see more below).

The results are processed into TSV and markdown reports. We use a custom Jekyll layout to present these reports as a static web-site with data elements made interactive through javascript. Tables are handled by datatables.js and plots are handled by plotly.js.

All code is available [on GitHub](https://github.com/daniel-munro/rat-twas-hub), originally forked from [gusevlab/TWAS_HUB](https://github.com/gusevlab/TWAS_HUB). All data for each trait is available from the <i class="far fa-file-archive" aria-hidden="true"></i> links in the [TRAITS](/traits){: .border} tab.

#### What does a TWAS association really mean?

Please read this [blog post](http://sashagusev.github.io/2017-10/twas-vulnerabilities.html) for much more about interpreting TWAS signals and the relationship between TWAS, other methods, and complex trait architectures.

#### Predictive models and weights

The predictive models for all analyses were generated using the default [Pantry (Pan-transcriptomic phenotyping) pipeline](https://github.com/PejLab/Pantry). Weights were generated from all available RNA phenotypes, and multiple RNA phenotype PCs and genotype PCs were used as covariates. See the [RatGTEx Portal](https://ratgtex.org) for information on the RNA-seq data sources and processing.

#### Interpretation of low significance models

All analyses include weakly predictive models up to a heritability P-value of 0.01. This means you will sometimes see models with negative cross-validation (adjusted) R<sup>2</sup> values because the heritable signal is not predictive after reducing 4/5 folds. These models are included primarily for individual gene look-up where the multiple-testing burden is negligible and weakly significant models may still be informative (alternatively, if you don't see a model for a gene it's because there wasn't a hint of signal in the data). For genomewide scans we recommend interpreting these models with caution.

#### Interpretation of the conditional analysis

The conditional analysis is a simple summary based step-wise model selection process that iteratively adds predictors to the model in decreasing order of conditional TWAS significance until no significant associations remain. Across models, conditional results should be interpreted as estimating the number of jointly significant models, but the selected models are not necessarily more likely to be causal than unselected features (either due to high correlation or different levels of noise). Rather, we recommend using a formal fine-mapping procedure (e.g. [FOCUS](https://github.com/mancusolab/ma-focus)).

Additionally, the SNP conditioning analysis (and Manhattan plots) provides an estimate of variance in the locus explained by the predicted model. A small fraction of variance explained is a strong indicator that the predicted model is tagging another causal feature (or there are multiple causal features in the locus). A large fraction of variance explained is consistent with the predicted model explaining all of the genetic effect - necessary but not sufficient for this to be the single causal mediator.

The conditional analysis uses an LD-reference panel and is therefore approximate, so you may see loci that behave unusually (for example, becoming extremely significant after conditioning). These are most likely instances of LD mismatch between reference and GWAS data.

#### Interpretation of coloc posteriors

All transcriptome-wide significant associations are run through the `coloc` colocalization model, with posterior probabilities PP3 (distinct causal variant) and PP4 (shared causal variant) reported in the locus view. `coloc` assumes a single causal variant model while TWAS directly models multiple xQTLs so we tend to use low PP3 as an indicator of colocalization rather than high PP4 (as done in [Raj et al. 2018 *Nature Genetics*](https://www.nature.com/articles/s41588-018-0238-1)).

# Acknowledgements

HS Rat RNA-seq datasets used for the transcriptomic models were produced by multiple studies coordinated by the [NIDA Center of Excellence for Genetics, Genomics, and Epigenetics of Substance Use Disorders in Outbred Rats](https://ratgenes.org). See the [RatGTEx homepage](https://ratgtex.org) for a list of investigators, funding, and data access. GWAS summary data was also processed by the center for the [projects listed here](/projects).

#### What else can I read about TWAS?

| [Munro et al. 2024 Nat Commun](https://pubmed.ncbi.nlm.nih.gov/39613793/) | Use of multiple RNA modalities |
| [Mancuso et al. 2018 Nat Genet](https://pubmed.ncbi.nlm.nih.gov/30926970/) | A method for fine-mapping credible sets of TWAS genes |
| [Barfield et al. 2018 Genet Epidemiol](https://pubmed.ncbi.nlm.nih.gov/29808603/) | A method for distinguishing co-localization in TWAS tests |
| [Gusev et al. 2018 biorxiv](https://doi.org/10.1101/330613) | TWAS of ovarian cancer |
| [Mancuso et al. 2018 Nat Commun](https://pubmed.ncbi.nlm.nih.gov/30287866/) | TWAS of prostate cancer |
| [Wu et al. 2018 Nat Genet](https://pubmed.ncbi.nlm.nih.gov/29915430/) | TWAS of breast cancer |
| [Gusev et al. 2018 Nat Genet](https://pubmed.ncbi.nlm.nih.gov/29632383/) | Integration of TWAS with chromatin features |
| [Mancuso et al. 2017 AJHG](https://pubmed.ncbi.nlm.nih.gov/28238358/) | TWAS of 30 traits and methods for cross-trait analyses |
| [Gusev et al. 2016 Nat Genet](https://pubmed.ncbi.nlm.nih.gov/26854917/) | Primary TWAS method paper |
{: .table}

## Change Log

| 2025-09-23 | Upgraded data to RatGTEx v3 (RefSeq annotations, new tissues) |
| 2025-09-08 | Added tag-based trait filtering and GWAS source publications |
| 2025-05-09 | Improved UI and table info, added cross-species gene search |
| 2025-02-12 | Added traits (281 total), improved trait metadata, and merged same-tissue datasets |
| 2024-09-03 | Added page with table of source projects |
| 2024-07-26 | Upgraded reference genome to mRatBN7.2 and added traits (123 total) for all RatGTEx tissues |
| 2024-05-31 | Added four RNA modalities that have multiple phenotypes per gene |
| 2024-04-25 | Initial beta release with 18 traits and expression and RNA stability from whole brain tissue |
{: .table}
