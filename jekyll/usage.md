---
title: Usage
permalink: usage/
---

# Usage

This page describes an example analysis of an individual trait. First, go to the [traits]({{ site.baseurl }}traits/) tab and search for "Locomotion velocity, NPP test", or follow this link to the ["Locomotion velocity, NPP test" trait view]({{ site.baseurl }}traits/novelty_seeking_test_total_velocity).

### Trait view

The first table of the Trait View ("Significant Loci") groups transcriptome-wide significant associations into contiguous loci after Bonferroni correction for all tested models. The rightmost column lists genes with at least one RNA model retained by stepwise conditional TWAS. The conditional GWAS P is evaluated at the original lead SNP, the SNP with the lowest marginal GWAS P in the locus. [Models tested]({{ site.baseurl }}models/)

**Shared xTWAS signals** lists other traits with at least one matching RNA model that passes the transcriptome-wide cutoff in both traits. For each trait, the cutoff is P < 0.05 divided by the total number of models in the release. A model matches when its weight file and full RNA phenotype ID are identical, so different tissues, modalities, or RNA phenotypes for the same gene do not count as a match. A gene is counted once even if several of its models match and pass. The table is sorted by shared gene count, then matching model count, and you can search by trait or gene. These counts use all models passing the cutoff before stepwise selection. Shared signals show association with both traits and leave the causal mechanism unresolved.

Expand **Directional secondary-trait test** to see a separate comparison. It starts with all transcriptome-wide significant models for the current trait and tests those same models in each secondary trait at P < 0.05 divided by the number of significant models in the current trait. Switching which trait is primary can change these results. The χ² ratio divides the mean squared Z-score of passing models with available secondary-trait Z-scores by the secondary trait's overall mean squared Z-score. The overall mean uses all rows in the secondary trait's results table, with missing Z-scores contributing zero. The correlation is Pearson's r across paired, signed model Z-scores that pass the directional test. The table shows a dash for r and its P-value when fewer than four finite Z-score pairs are available or when either set of Z-scores is constant.

The "Associations by Panel" table shows the breakdown of associations by transcriptomic panel. These are ordered by the average TWAS chi<sup>2</sup> statistic in the panel -- an estimate of the average trait heritability explained by predictors from that panel. The columns also report the # and % of significant associations from that study.

### Locus view

Open locus 2 in the "Locomotion velocity, NPP test" associations table to see the GWAS results before and after conditioning on the selected RNA models. The summary above the plot reports the marginal and conditional GWAS P at the same original lead SNP. [Locus 2]({{ site.baseurl }}traits/novelty_seeking_test_total_velocity/2/)

The **Associated models** table reports two separate results for follow-up. **Joint selection** records whether stepwise conditional TWAS retained the RNA model alongside the other selected models at the locus. **Shared-variant posterior (PP4)** measures support for a single causal variant shared by the RNA phenotype and the GWAS trait under the coloc model. **Distinct-variant posterior (PP3)** measures support for each having its own causal variant. A causal role for the RNA phenotype requires additional evidence. [About the analyses]({{ site.baseurl }}about)

### Gene view

Click on [Ctsc]({{ site.baseurl }}genes/Ctsc/) to go to the Gene View. The top table ("Models") shows all of the predictive models that have been computed for this gene and their respective performance.

The second table of the gene view ("Trait Associations") shows a heatmap of association for this gene between all traits (rows) and all models. We order the heatmap by the "Avg chi<sup>2</sup> ratio" column, which is computed as the average chi<sup>2</sup> for the gene-trait pair (across all models) divided by the average chi<sup>2</sup> for all genes in the listed trait (across all models). This normalization accounts for sample size and heritability differences between traits and emphasize associations that are stronger than expected by chance (without the normalization, highly heritable and polygenic traits like height, for example, would constantly be at the top of the list simply because they have so many detectable causal variants). The subsequent columns list the raw average chi<sup>2</sup> statistic, maximum chi<sup>2</sup> statistic across all models (to filter for model-specific associations), and then the individual Z-scores for each model.
