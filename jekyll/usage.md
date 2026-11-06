---
title: Usage
permalink: usage/
---

# Usage

This page describes an example analysis of an individual trait. First, go to the [traits]({{ site.baseurl }}traits/) tab and search for "Locomotion velocity, NPP test", or follow this link to the ["Locomotion velocity, NPP test" trait view]({{ site.baseurl }}traits/novelty_seeking_test_total_velocity).

### Trait view

The trait header lists the GWAS sample size and source study and links to the full xTWAS results. Use **Significant xTWAS models** to filter the transcriptome-wide hits by gene, tissue, RNA modality, or locus. Each row represents one RNA model. The **Significant loci** table groups hits into contiguous regions, with links to each locus and to the genes with an RNA model retained by stepwise conditional TWAS. The conditional GWAS P is evaluated at the original lead SNP, the SNP with the lowest marginal GWAS P in the locus. [Models tested]({{ site.baseurl }}models/)

**Shared xTWAS signals** lists other traits with at least one matching RNA model that passes the transcriptome-wide cutoff in both traits. For each trait, the cutoff is P < 0.05 divided by the total number of models in the release. A model matches when its weight file and full RNA phenotype ID are identical, so different tissues, modalities, or RNA phenotypes for the same gene do not count as a match. A gene is counted once even if several of its models match and pass. The table is sorted by shared gene count, then matching model count, and you can search by trait or gene. These counts use all models passing the cutoff before stepwise selection. Shared signals show association with both traits and leave the causal mechanism unresolved.

Expand **Directional secondary-trait test** to see a separate comparison. It starts with all transcriptome-wide significant models for the current trait and tests those same models in each secondary trait at P < 0.05 divided by the number of significant models in the current trait. Switching which trait is primary can change these results. The χ² ratio divides the mean squared Z-score of passing models with available secondary-trait Z-scores by the secondary trait's overall mean squared Z-score. The overall mean uses all rows in the secondary trait's results table, with missing Z-scores contributing zero. The correlation is Pearson's r across paired, signed model Z-scores that pass the directional test. The table shows a dash for r and its P-value when fewer than four finite Z-score pairs are available or when either set of Z-scores is constant.

**Associations by panel** counts distinct genes with hits and significant models among models with a TWAS P-value for each tissue and modality.

### Locus view

Open locus 2 in the "Locomotion velocity, NPP test" associations table to compare marginal and conditional GWAS results on the same genomic and −log<sub>10</sub> P scales. The summary above the plots reports the marginal and conditional GWAS P at the same original lead SNP. [Locus 2]({{ site.baseurl }}traits/novelty_seeking_test_total_velocity/2/)

**Associated models** groups RNA features by gene and shows tissue, modality, signed TWAS Z, TWAS P, and cross-validation R². **Joint selection** records whether stepwise conditional TWAS retained the RNA model. The colocalization posteriors provide separate evidence for follow-up. **PP4** measures support for a shared causal variant and **PP3** measures support for distinct causal variants under the coloc model. Expand a model's diagnostics to see its full RNA phenotype ID, which is also included in the locus TSV download. [About the analyses]({{ site.baseurl }}about)

### Gene view

Open the gene page to see the strongest significant RNA model for each trait, ranked by TWAS P. The lead table includes its tissue, modality, and a link to the locus when available. The RatGTEx link opens expression and xQTL results. The **RNA model summary** shows the selected predictor and cross-validation R² recorded in the reference TWAS results. [Ctsc gene page]({{ site.baseurl }}genes/Ctsc/)

Expand **Model diagnostics by method** to compare heritability and predictor performance. Expand **Full trait association heatmap and χ² rankings** to see model-level TWAS Z-scores across traits and the average and maximum χ² statistics. The average χ² ratio compares the gene's mean squared Z-score with the overall mean squared Z-score for that trait.
