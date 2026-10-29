# Compute directed trait comparisons while parsing each association table twice.
suppressPackageStartupMessages(library(tidyverse))

pheno_to_gene_id <- function(pheno_id) str_replace(pheno_id, "__.+$", "")
tbl_traits <- read_tsv("data/traits.tsv", col_types = "ccicccc") |>
  left_join(read_tsv("data/traits.summary.tsv", col_types = "ciiid"),
            by = "ID", relationship = "one-to-one")
n_models <- length(read_lines("data/all_models.tsv")) - 1

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
  read_assocs(path) |> filter(TWAS.P < 0.05 / n_models)
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
  top_models <- seq_len(nrow(trait_assocs))
  fout_pleio <- str_glue("jekyll/_data/trait_pleio/{trait}.tsv")
  # ---- Get pleiotropic associations
  df_pleiot <- tibble(
    trait_id = character(),
    chisq_ratio = numeric(),
    num_genes = integer(),
    num_genes_twas = integer(),
    pct_genes_twas = numeric(),
    corr = numeric(),
    p_val = numeric(),
    genes = character()
  )

  for (i_trait2 in seq_len(nrow(tbl_traits))) {
    if ((i_trait2 != i_trait) && length(top_models) > 0) {
      trait2_assocs <- secondary_assocs[[i_trait2]]
      pleio <- inner_join(trait_assocs, trait2_assocs, by = "model_index", relationship = "one-to-one") |>
        filter(TWAS.P.x < 0.05 / n_models,
               TWAS.P.y < 0.05 / length(top_models))
      if (nrow(pleio) > 0) {
        if(nrow(pleio) >= 4) {
          tst <- with(pleio, cor.test(TWAS.Z.x, TWAS.Z.y))
        } else {
          tst <- tibble(est = 0, p.value = 1)
        }
        genes <- sort(unique(pheno_to_gene_id(pleio$ID)))
        n_genes_twas <- pleio |>
          filter(TWAS.P.y < 0.05 / n_models) |>
          with(length(unique(pheno_to_gene_id(ID))))
        df_pleiot <- bind_rows(df_pleiot, tibble(
          trait_id = tbl_traits$ID[i_trait2],
          chisq_ratio = mean(pleio$TWAS.Z.y^2, na.rm = TRUE) / tbl_traits$AVG.CHISQ[i_trait2],
          num_genes = length(genes),
          num_genes_twas = n_genes_twas,
          pct_genes_twas = round(100 * n_genes_twas / tbl_traits$NUM.JOINT.GENES[i_trait], 1),
          corr = tst$est,
          p_val = tst$p.value,
          genes = str_c(genes, collapse = ",")
        ))
      }
    }
  }

  df_pleiot |>
    replace_na(list(pct_genes_twas = 0)) |>
    select(
      trait_id,
      chisq_ratio,
      num_genes,
      num_genes_twas,
      pct_genes_twas,
      corr,
      p_val,
      genes,
    ) |>
    mutate(
      chisq_ratio = round(chisq_ratio, 2),
      corr = round(corr, 2),
      p_val = format(p_val, digits = 3, scientific = TRUE),
    ) |>
    write_tsv(fout_pleio)
}
