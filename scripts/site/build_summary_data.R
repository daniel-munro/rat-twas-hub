# Load the data and generate summary stats and other reference data.

# Inputs:
# data/traits.par
# data/traits.par.nfo
# data/all_models.par
# data/genes_n_assoc.nfo
# data/genes_n_models.nfo
# data/panels.par

# Outputs:
# jekyll/_data/traits.tsv
# jekyll/_data/genes.json
# jekyll/_data/stats.yml
# jekyll/_data/panels.tsv

suppressPackageStartupMessages(library(tidyverse))
library(yaml)

gene_id <- function(pheno_id) {
  # Extract gene ID from RNA phenotype ID
  stringr::str_replace(pheno_id, "__.+$", "")
}

## Print trait index

traits_nfo <- read_tsv("data/traits.par.nfo", col_types = "ciiid")

tbl_traits <- read_tsv("data/traits.par", col_types = "ccicccc") |>
  left_join(traits_nfo, by = "ID")

tbl_traits |>
  select(ID, N, PROJECT, TAGS, NAME, DESCRIPTION, NUM_LOCI = NUM.LOCI, NUM_JOINT_GENES = NUM.JOINT.GENES, NUM_GENES = NUM.GENES) |>
  write_tsv("jekyll/_data/traits.tsv")

## Make genes.json

all_genes <- read_tsv("data/all_models.par", col_types = cols(ID = "c", .default = "-")) |>
  mutate(gene_id = gene_id(ID)) |>
  with(sort(unique(gene_id)))
genes_n_assoc <- read_delim("data/genes_n_assoc.nfo", delim = " ", col_types = "ci", col_names = FALSE) |>
  deframe()
genes_n_models <- read_delim("data/genes_n_models.nfo", delim = " ", col_types = "ci", col_names = FALSE) |>
  deframe()

df_genes <- tibble(gene_id = all_genes) |>
  mutate(
    n_assoc = genes_n_assoc[gene_id],
    n_models = genes_n_models[gene_id],
    gene_link = str_glue('<em><a href=\\"./{gene_id}\\">{gene_id}</a></em>'),
  ) |>
  replace_na(list(n_assoc = 0, n_models = 0))

cat('{\n"data":[\n',
    str_c(str_glue('["{df_genes$gene_link}",{df_genes$n_assoc},{df_genes$n_models}]'), collapse = ",\n"),
    "]\n}",
    sep = "",
    file = "jekyll/genes.json")

## Stats for trait and gene index pages

gene_stats <- list(
  n_traits = nrow(tbl_traits),
  n_loci = formatC(sum(tbl_traits$NUM.LOCI), format = "f", big.mark = ",", drop0trailing = TRUE),
  n_gene_trait_assocs = formatC(sum(tbl_traits$NUM.GENES), format = "f", big.mark = ",", drop0trailing = TRUE),
  n_genes = formatC(nrow(df_genes), format = "f", big.mark = ",", drop0trailing = TRUE),
  n_models = formatC(sum(df_genes$n_models), format = "f", big.mark = ",", drop0trailing = TRUE)
)
dir.create("jekyll/_data", showWarnings = FALSE, recursive = TRUE)
yaml::write_yaml(gene_stats, file = "jekyll/_data/stats.yml")

## Make panels.tsv with model counts

panels <- readr::read_tsv("data/panels.par", col_types = "cccccc")

panel_counts <- read_tsv(
  "data/all_models.par", col_types = cols(PANEL = "c", .default = "-")
) |>
  count(PANEL, name = "N_MODELS")

panels_out <- panels |>
  left_join(panel_counts, by = "PANEL", relationship = "one-to-one")

write_tsv(panels_out, "jekyll/_data/panels.tsv")
