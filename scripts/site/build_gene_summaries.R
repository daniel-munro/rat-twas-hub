# Generate gene leads separately from annotations and the full TWAS heatmap.
suppressPackageStartupMessages(library(tidyverse))

genes <- commandArgs(trailingOnly = TRUE)
traits <- read_tsv("data/traits.tsv", col_types = "ccicccc")
panels <- read_tsv("data/panels.tsv", col_types = "-c-cc-")
models <- read_tsv("data/all_models.tsv",
  col_types = cols(PANEL = "c", ID = "c", .default = "-")) |>
  mutate(gene_id = str_replace(ID, "__.+$", ""))
if (length(genes) > 0) models <- filter(models, gene_id %in% genes)

# Predictor selection uses cross-validation of the RNA models and is shared
# across traits. Use the recorded choice and R² from the reference TWAS table.
predictors <- read_tsv(traits$OUTPUT[[1]],
  col_types = cols(PANEL = "c", ID = "c", MODEL = "c", MODELCV.R2 = "c", .default = "-")) |>
  filter(!is.na(MODEL), MODEL != "") |>
  distinct(PANEL, ID, MODEL, MODELCV.R2)
models <- models |>
  left_join(predictors, by = c("PANEL", "ID"), relationship = "one-to-one") |>
  left_join(panels, by = "PANEL", relationship = "many-to-one") |>
  arrange(TISSUE, MODALITY, ID) |>
  select(gene_id, tissue = TISSUE, modality = MODALITY, pheno_id = ID,
         predictor = MODEL, cv_r2 = MODELCV.R2)

hits <- map2(traits$ID, traits$NAME, function(trait, name) {
  read_tsv(str_glue("jekyll/_data/trait_hits/{trait}.tsv"),
    col_types = "ccccdcc") |>
    mutate(trait_id = trait, trait_name = name, twas_p = as.numeric(twas_p))
}) |>
  list_rbind() |>
  semi_join(distinct(models, gene_id), by = "gene_id") |>
  arrange(twas_p, desc(abs(twas_z)), trait_id, tissue, modality, pheno_id, locus) |>
  distinct(gene_id, trait_id, .keep_all = TRUE) |>
  mutate(twas_p = format(twas_p, digits = 3, scientific = TRUE, trim = TRUE),
         twas_z = round(twas_z, 2))

hits_by_gene <- split(select(hits, -gene_id), hits$gene_id)
models_by_gene <- split(select(models, -gene_id), models$gene_id)
output <- "jekyll/_data/gene_summaries"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
for (gene in names(models_by_gene)) {
  gene_hits <- hits_by_gene[[gene]]
  jsonlite::write_json(
    list(hits = if (is.null(gene_hits)) list() else gene_hits, models = models_by_gene[[gene]]),
    file.path(output, paste0(gene, ".json")), na = "null"
  )
}
