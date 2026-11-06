# Compute shared xTWAS signals and directional tests in two passes over the inputs.
suppressPackageStartupMessages(library(tidyverse))

pheno_to_gene_id <- function(pheno_id) str_replace(pheno_id, "__.+$", "")
tbl_traits <- read_tsv("data/traits.tsv", col_types = "ccicccc") |>
  left_join(read_tsv("data/traits.summary.tsv", col_types = "ciiid"),
            by = "ID", relationship = "one-to-one")
n_models <- length(read_lines("data/all_models.tsv")) - 1
twas_threshold <- 0.05 / n_models

read_assocs <- function(path) {
  df <- read_tsv(path, num_threads = 1,
    col_types = cols(FILE = "c", ID = "c", TWAS.Z = "d", TWAS.P = "d", .default = "-")) |>
    filter(!is.na(TWAS.P))
  # The original one-to-one joins require unique model keys for each trait.
  stopifnot(!anyDuplicated(df[c("FILE", "ID")]))
  df
}

# First pass: significant primary models and their union across configured traits.
top_assocs <- lapply(tbl_traits$OUTPUT, function(path) {
  cat("Selecting primary models:", path, "\n")
  read_assocs(path) |> filter(TWAS.P < twas_threshold)
})
model_keys <- bind_rows(top_assocs) |>
  distinct(FILE, ID) |>
  mutate(model_index = row_number())
top_assocs <- lapply(top_assocs, function(df) {
  inner_join(df, model_keys, by = c("FILE", "ID"), relationship = "one-to-one") |>
    select(model_index, ID, TWAS.Z, TWAS.P)
})

# Second pass: keep only these union models in each secondary trait. Integer keys
# avoid retaining repeated file paths and phenotype IDs across the trait tables.
secondary_assocs <- lapply(tbl_traits$OUTPUT, function(path) {
  cat("Reading secondary models:", path, "\n")
  inner_join(read_assocs(path), model_keys, by = c("FILE", "ID"),
             relationship = "one-to-one") |>
    select(model_index, TWAS.Z, TWAS.P)
})
cat("Retained", nrow(model_keys), "union models for", nrow(tbl_traits), "traits.\n")

for (i_trait in seq_len(nrow(tbl_traits))) {
  trait <- tbl_traits$ID[i_trait]
  cat("Comparing:", trait, "\n")
  trait_assocs <- top_assocs[[i_trait]]
  n_primary_hits <- nrow(trait_assocs)
  comparisons <- tibble(
    trait_id = character(),
    num_genes_twas = integer(),
    num_models_twas = integer(),
    genes_twas = character(),
    chisq_ratio = numeric(),
    num_genes = integer(),
    corr = numeric(),
    p_val = numeric(),
    genes = character()
  )

  for (i_trait2 in seq_len(nrow(tbl_traits))) {
    if (i_trait2 == i_trait || n_primary_hits == 0) next

    matched <- inner_join(trait_assocs, secondary_assocs[[i_trait2]],
                          by = "model_index", relationship = "one-to-one")
    shared <- filter(matched, TWAS.P.y < twas_threshold)
    directional <- filter(matched, TWAS.P.y < 0.05 / n_primary_hits)
    if (nrow(shared) == 0 && nrow(directional) == 0) next

    shared_genes <- sort(unique(pheno_to_gene_id(shared$ID)))
    directional_genes <- sort(unique(pheno_to_gene_id(directional$ID)))
    z_pairs <- filter(directional, is.finite(TWAS.Z.x), is.finite(TWAS.Z.y))
    corr <- p_val <- NA_real_
    if (nrow(z_pairs) >= 4 &&
        n_distinct(z_pairs$TWAS.Z.x) > 1 && n_distinct(z_pairs$TWAS.Z.y) > 1) {
      tst <- with(z_pairs, cor.test(TWAS.Z.x, TWAS.Z.y))
      corr <- unname(tst$estimate)
      p_val <- tst$p.value
    }
    comparisons <- bind_rows(comparisons, tibble(
      trait_id = tbl_traits$ID[i_trait2],
      num_genes_twas = length(shared_genes),
      num_models_twas = nrow(shared),
      genes_twas = paste(shared_genes, collapse = ","),
      chisq_ratio = mean(directional$TWAS.Z.y^2, na.rm = TRUE) / tbl_traits$AVG.CHISQ[i_trait2],
      num_genes = length(directional_genes),
      corr = corr,
      p_val = p_val,
      genes = paste(directional_genes, collapse = ",")
    ))
  }

  comparisons |>
    mutate(
      chisq_ratio = round(chisq_ratio, 2),
      corr = round(corr, 2),
      p_val = signif(p_val, 3),
    ) |>
    write_tsv(str_glue("jekyll/_data/trait_pleio/{trait}.tsv"), na = "")
}
