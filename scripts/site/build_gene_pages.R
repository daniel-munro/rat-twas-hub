# Generate gene pages

# Inputs:
# data/panels.par
# data/all_models.par
# data/traits.par
# data/traits.par.nfo
# data/twas_out/{trait}.all.tsv (for each trait listed in traits.par)

# Outputs:
# jekyll/genes/{gene}.md
# jekyll/_data/gene_models/{gene}.tsv

suppressPackageStartupMessages(library(tidyverse))

gene_id <- function(pheno_id) {
  # Extract gene ID from RNA phenotype ID
  stringr::str_replace(pheno_id, "__.+$", "")
}

shorten_ids <- function(pheno_id, modality) {
  # Remove extra information from RNA phenotype IDs
  # alternative TSS/polyA:            {gene_id}__grp_{grp#}_{up|down}stream_{transcript_id}      -> {transcript_id}
  # gene expression & mRNA stability: {gene_id}                                                  -> {gene_id}
  # intron excision ratio:            {gene_id}__{chrom}_{pos1}_{pos2}_clu_{cluster_id}_{strand} -> {chrom}_{pos1}_{pos2}
  # isoform ratio:                    {gene_id}__{transcript_id}                                 -> {transcript_id}
  dplyr::case_when(
    modality %in% c("alternative TSS", "alternative polyA") ~ pheno_id |> stringr::str_replace("^.+stream_", ""),
    modality == "intron excision ratio" ~ pheno_id |> stringr::str_replace("^.+__", "") |> stringr::str_replace("_clu_.+$", ""),
    modality == "isoform ratio" ~ pheno_id |> stringr::str_replace("^.+__", ""),
    .default = pheno_id
  )
}

traits_nfo <- read_tsv("data/traits.par.nfo", col_types = "ciiid")

tbl_traits <- read_tsv("data/traits.par", col_types = "ccicccc") |>
  mutate(link = str_glue("<a href='{{{{ site.baseurl }}}}traits/{ID}'>{NAME}</a>")) |>
  left_join(traits_nfo, by = "ID", relationship = "one-to-one")

cat("Reading all trait TWAS outputs\n")

tbl_zscores <- tbl_traits |>
  rename(trait_id = ID) |>
  reframe(
    read_tsv(OUTPUT, col_types = cols(PANEL = "c", ID = "c", TWAS.Z = "d", TWAS.P = "d", .default = "-")) |>
      rename(pheno_id = ID, z = TWAS.Z, p = TWAS.P),
    .by = trait_id
  ) |>
  mutate(gene_id = gene_id(pheno_id)) |>
  replace_na(list(z = 0, p = 1))

model_order <- tbl_zscores |>
  summarise(
    sum_z_sq = sum(z^2, na.rm = TRUE),
    .by = c(gene_id, PANEL, pheno_id)
  ) |>
  arrange(gene_id, desc(sum_z_sq)) |>
  mutate(model_num = seq_len(n()), .by = gene_id)

cat(str_glue("Loaded {scales::comma(nrow(tbl_zscores))} trait-model rows for {length(unique(tbl_zscores$trait_id))} traits."),"\n")

tbl_panels <- read_tsv("data/panels.par", col_types = "-c-cc-")

tbl_models <- read_tsv("data/all_models.par", col_types = "cccciiiddddddddddddd") |>
  rename(pheno_id = ID) |>
  left_join(tbl_panels, by = "PANEL", relationship = "many-to-one") |>
  mutate(gene_id = gene_id(pheno_id)) |>
  left_join(model_order, by = c("gene_id", "PANEL", "pheno_id"), relationship = "one-to-one") |>
  arrange(gene_id, model_num) |>
  select(
    gene_id,
    model_num,
    tissue = TISSUE,
    modality = MODALITY,
    pheno_id,
    hsq = HSQ,
    hsq_se = HSQ.SE,
    hsq_pv = HSQ.PV,
    eqtl_r2 = TOP1.R2,
    blup_r2 = BLUP.R2,
    enet_r2 = ENET.R2,
    lasso_r2 = LASSO.R2,
    eqtl_pv = TOP1.PV,
    blup_pv = BLUP.PV,
    enet_pv = ENET.PV,
    lasso_pv = LASSO.PV,
  ) |>
  mutate(
    pheno_id = shorten_ids(pheno_id, modality),
    hsq = round(hsq, 2),
    hsq_se = round(hsq_se, 2),
    hsq_pv = format(hsq_pv, digits = 3, scientific = TRUE, trim = TRUE),
    eqtl_r2 = round(eqtl_r2, 2),
    blup_r2 = round(blup_r2, 2),
    enet_r2 = round(enet_r2, 2),
    lasso_r2 = round(lasso_r2, 2),
    eqtl_pv = format(eqtl_pv, digits = 3, scientific = TRUE, trim = TRUE),
    blup_pv = format(blup_pv, digits = 3, scientific = TRUE, trim = TRUE),
    enet_pv = format(enet_pv, digits = 3, scientific = TRUE, trim = TRUE),
    lasso_pv = format(lasso_pv, digits = 3, scientific = TRUE, trim = TRUE),
  )

tbl_genes_traits <- tbl_zscores |>
  left_join(
    select(model_order, gene_id, PANEL, pheno_id, model_num),
    by = c("gene_id", "PANEL", "pheno_id"),
    relationship = "many-to-one"
  ) |>
  arrange(gene_id, trait_id, model_num) |>
  mutate(
    significant = p < 0.05 / nrow(tbl_models),
    z_html = if_else(
      significant,
      sprintf("<td><b>%2.1f</b></td>", round(z, 1)),
      sprintf("<td>%2.1f</td>", round(z, 1))
    )
  ) |>
  summarise(
    gene_avg_chisq = mean(z^2, na.rm = TRUE),
    gene_max_chisq = max(z^2, na.rm = TRUE),
    z_cells = paste(z_html, collapse = ""),
    .by = c(gene_id, trait_id)
  ) |>
  left_join(
    tbl_traits |>
      select(trait_id = ID, trait_link = link, trait_avg_chisq = AVG.CHISQ),
    by = "trait_id",
    relationship = "many-to-one"
  ) |>
  mutate(
    avg_chisq_ratio = if_else(trait_avg_chisq != 0, gene_avg_chisq / trait_avg_chisq, NA_real_),
    row_html = sprintf(
      "<tr><td>%s</td><td>%2.1f</td><td>%2.1f</td><td>%2.1f</td>%s</tr>",
      trait_link, round(avg_chisq_ratio, 1), round(gene_avg_chisq, 1), round(gene_max_chisq, 1), z_cells
    )
  )

all_genes <- sort(unique(tbl_models$gene_id))

for (i in seq_along(all_genes)) {
  g_id <- all_genes[i]
  fout_page <- str_glue("jekyll/genes/{g_id}.md")
  fout_models <- str_glue("jekyll/_data/gene_models/{g_id}.tsv")

  system(str_glue("mkdir -p jekyll/genes"))
  system(str_glue("mkdir -p jekyll/_data/gene_models"))

  c(
    "---",
    str_glue('title: "{g_id}"'),
    str_glue("permalink: genes/{g_id}/"),
    "layout: gene",
    str_glue("id: {g_id}"),
    "---"
  ) |>
    write_lines(fout_page)

  # Trait association rows (columns ordered by the same precomputed model order)
  tbl_genes_traits |>
    filter(gene_id == g_id) |>
    pull(row_html) |>
    write(file = fout_page, append = TRUE)

  tbl_models |>
    filter(gene_id == g_id) |>
    select(-gene_id) |>
    write_tsv(fout_models)

  if (i %% 1000 == 0) cat(i, "genes processed\n")
}

cat("Finished writing gene pages with trait associations.\n")
