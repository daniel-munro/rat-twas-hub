# Generate gene pages

# Inputs:
# data/gene_names.tsv
# data/panels.par
# data/all_models.par
# data/traits.par
# data/traits.par.nfo

# Outputs:
# jekyll/genes/{gene}.md
# jekyll/_data/gene_models/{gene}.tsv

suppressPackageStartupMessages(library(tidyverse))

gene_id <- function(pheno_ID) {
  # Extract gene ID from RNA phenotype ID
  str_extract(pheno_ID, "^[^:.]+")
}

shorten_ids <- function(ID, modality) {
  # Remove extra information from RNA phenotype IDs
  case_when(
    modality %in% c("alternative polyA", "alternative TSS", "isoform ratio") ~ str_extract(ID, "[^:.]+$"),
    modality == "intron excision ratio" ~ str_c("chr", ID |> str_replace("^[^:]+:", "") |> str_replace(":[^:]+$", "")),
    .default = ID
  )
}

tbl_panels <- read_tsv("data/panels.par", col_types = "ccccci")

tbl_models <- read_tsv("data/all_models.par", col_types = "ccciiiiddddddddddddd") |>
  left_join(select(tbl_panels, PANEL, TISSUE, MODALITY), by = "PANEL") |>
  mutate(gene_id = gene_id(ID),
         ID = shorten_ids(ID, MODALITY))

traits_nfo <- read_tsv("data/traits.par.nfo", col_types = "ciiid")

tbl_traits <- read_tsv("data/traits.par", col_types = "ccicccc") |>
  mutate(link = str_glue("<a href='{{{{ site.baseurl }}}}traits/{ID}'>{NAME}</a>")) |>
  left_join(traits_nfo, by = "ID")

# iterate and make each gene file
all_genes <- sort(unique(tbl_models$gene_id))

# Load gene names since Ensembl IDs were used for TWAS
gene_names <- read_tsv("data/gene_names.tsv", col_types = "cc") |>
  complete(ID = all_genes) |>
  mutate(NAME = if_else(is.na(NAME), ID, NAME)) |>
  deframe()

for (i in seq_along(all_genes)) {
  g_id <- all_genes[i]
  fout_page <- str_glue("jekyll/genes/{g_id}.md")
  fout_models <- str_glue("jekyll/_data/gene_models/{g_id}.tsv")
  
  system(str_glue("mkdir -p jekyll/_data/gene_models"))
  
  gene_models <- tbl_models |>
    filter(gene_id == g_id) |>
    mutate(num = seq_len(n()))
  
  gene_models |>
    select(
      num,
      tissue = TISSUE,
      modality = MODALITY,
      id = ID,
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
      hsq = round(hsq, 2),
      hsq_se = round(hsq_se, 2),
      hsq_pv = format(hsq_pv, digits = 3, scientific = TRUE),
      eqtl_r2 = round(eqtl_r2, 2),
      blup_r2 = round(blup_r2, 2),
      enet_r2 = round(enet_r2, 2),
      lasso_r2 = round(lasso_r2, 2),
      eqtl_pv = format(eqtl_pv, digits = 3, scientific = TRUE),
      blup_pv = format(blup_pv, digits = 3, scientific = TRUE),
      enet_pv = format(enet_pv, digits = 3, scientific = TRUE),
      lasso_pv = format(lasso_pv, digits = 3, scientific = TRUE),
    ) |>
    write_tsv(fout_models)
  
  c(
    "---",
    str_glue('title: "{gene_names[g_id]}"'),
    str_glue("permalink: genes/{g_id}/"),
    "layout: gene",
    str_glue("id: {g_id}"),
    "---"
  ) |>
    write_lines(fout_page)
  
  if (i %% 1000 == 0) cat(i, "\n")
}

# Iterate over each trait and count the number of significant genes
cat("Adding trait associations to gene pages\n")
for (i in seq_len(nrow(tbl_traits))) {
  system(str_glue('bash scripts/site/add_trait_to_gene_pages.sh {tbl_traits$OUTPUT[i]} "{tbl_traits$link[i]}" {tbl_traits$AVG.CHISQ[i]}\n'))
  cat(i, "\n")
}
