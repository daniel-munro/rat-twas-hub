---
title: Usage
permalink: usage/
---

# Usage

This page describes an example analysis of an individual trait. First, go to the [traits]({{ site.baseurl }}traits/) tab and search for "Locomotion velocity, NPP test", or follow this link to the ["Locomotion velocity, NPP test" trait view]({{ site.baseurl }}traits/novelty_seeking_test_total_velocity).

### Trait view

The first table of the Trait View ("Significant Loci") shows all of the transcriptome-wide significant associations for the given trait (after bonferroni correction for all [models]({{ site.baseurl }}models/) tested). Loci have been grouped into contiguous blocks and model selection run on each locus to identify the independently significant genes (which are reported in the right-most column).

The second table of the Trait View ("Pleiotropic Associations") shows all pleiotropic associations to other traits for any of the independently significant genes. You can filter this table by trait or gene. The table is ordered by the "chi<sup>2</sup> ratio" which is computed as the average chi<sup>2</sup> statistic for the selected genes in the secondary trait, divided by the average statistic for all genes in the secondary trait. Ignoring issues of LD, this is an estimate of the heritability enrichment of the target genes relative to all genes and tends to provide reasonable results.
The remaining columns list the number of significant genes in the target trait at Bonferroni correction [+] and at transcriptome-wide significance [++], the correlation of effect-sizes across the [+] genes, as well as links to each of the [+] genes.

_NB: (1) although we represent the associations as genes, we are actually testing for the same associated tissue/model; (2) we use the jointly significant genes so that the pleiotropy statistics are computed from independent data points, but other correlated genes from each locus should also be treated as pleiotropic_.

The third table of the Trait View ("Associations by Panel") shows the breakdown of associations by transcriptomic panel. These are ordered by the average TWAS chi<sup>2</sup> statistic in the panel -- an estimate of the average trait heritability explained by predictors from that panel. The columns also report the # and % of significant associations from that study.

### Locus view

Click on [locus #2]({{ site.baseurl }}traits/novelty_seeking_test_total_velocity/2/) in the "Locomotion velocity, NPP test" associations table to go to the Locus View for the _Ctsc_/_Grm5_/_Tyr_ locus. The top panel shows a Manhattan plot of the GWAS association before and after conditioning on the predicted transcriptomic variation (see [About]({{ site.baseurl }}about) for more details on conditioning). The next panel is a table ("Associated Models") that shows all of the significantly associated models, their model performance, and posterior probability of a single shared causal variant.

### Gene view

Click on [Ctsc]({{ site.baseurl }}genes/Ctsc/) to go to the Gene View. The top table ("Models") shows all of the predictive models that have been computed for this gene and their respective performance.

The second table of the gene view ("Trait Associations") shows a heatmap of association for this gene between all traits (rows) and all models. We order the heatmap by the "Avg chi<sup>2</sup> ratio" column, which is computed as the average chi<sup>2</sup> for the gene-trait pair (across all models) divided by the average chi<sup>2</sup> for all genes in the listed trait (across all models). This normalization accounts for sample size and heritability differences between traits and emphasize associations that are stronger than expected by chance (without the normalization, highly heritable and polygenic traits like height, for example, would constantly be at the top of the list simply because they have so many detectable causal variants). The subsequent columns list the raw average chi<sup>2</sup> statistic, maximum chi<sup>2</sup> statistic across all models (to filter for model-specific associations), and then the individual Z-scores for each model.
