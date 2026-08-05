# Load the data and generate summary stats and other reference data.

# Inputs:
# data/traits.tsv
# data/traits.summary.tsv
# data/all_models.tsv
# data/genes_n_assoc.tsv
# data/genes_n_models.tsv
# data/panels.tsv
# GTF annotation passed as the first command-line argument

# Outputs:
# jekyll/_data/traits.tsv
# jekyll/_data/genes.json
# jekyll/_data/stats.yml
# jekyll/_data/panels.tsv

suppressPackageStartupMessages(library(tidyverse))
library(yaml)

args <- commandArgs(trailingOnly = TRUE)
gtf_file <- args[[1]]

gtf_attributes <- function(attributes, key) {
  pattern <- as.character(str_glue("(?:^|; )\\Q{key}\\E \\\"([^\\\"]*)\\\""))
  stringr::str_match(attributes, stringr::regex(pattern))[, 2]
}

html_escape <- function(value) {
  str_replace_all(value, c("&" = "&amp;", "<" = "&lt;", ">" = "&gt;", "\\\"" = "&quot;", "'" = "&#39;"))
}

gene_id <- function(pheno_id) {
  # Extract gene ID from RNA phenotype ID
  stringr::str_replace(pheno_id, "__.+$", "")
}

## Print trait index

traits_nfo <- read_tsv("data/traits.summary.tsv", col_types = "ciiid")

tbl_traits <- read_tsv("data/traits.tsv", col_types = "ccicccc") |>
  left_join(traits_nfo, by = "ID")

tbl_traits |>
  select(ID, N, PROJECT, TAGS, NAME, DESCRIPTION, NUM_LOCI = NUM.LOCI, NUM_JOINT_GENES = NUM.JOINT.GENES, NUM_GENES = NUM.GENES) |>
  write_tsv("jekyll/_data/traits.tsv")

## Make genes.json

all_genes <- read_tsv("data/all_models.tsv", col_types = cols(ID = "c", .default = "-")) |>
  mutate(gene_id = gene_id(ID)) |>
  with(sort(unique(gene_id)))
genes_n_assoc <- read_tsv("data/genes_n_assoc.tsv", col_types = "ci", col_names = FALSE) |>
  deframe()
genes_n_models <- read_tsv("data/genes_n_models.tsv", col_types = "ci", col_names = FALSE) |>
  deframe()

gene_names <- read_tsv(
  gtf_file,
  comment = "#",
  col_names = c("chrom", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"),
  col_types = "--c-----c"
) |>
  filter(feature == "gene") |>
  transmute(
    gene_id = gtf_attributes(attributes, "gene_id"),
    full_name = gtf_attributes(attributes, "description")
  )

df_genes <- tibble(gene_id = all_genes) |>
  left_join(gene_names, by = "gene_id", relationship = "one-to-one") |>
  mutate(
    n_assoc = genes_n_assoc[gene_id],
    n_models = genes_n_models[gene_id],
    gene_link = str_glue('<em><a href="./{gene_id}">{gene_id}</a></em>'),
    full_name_html = if_else(
      is.na(full_name),
      "",
      str_glue('<span class="gene-full-name" title="{html_escape(full_name)}">{html_escape(full_name)}</span>')
    )
  ) |>
  replace_na(list(n_assoc = 0, n_models = 0))

gene_rows <- pmap(
  select(df_genes, gene_link, full_name_html, n_assoc, n_models),
  \(gene_link, full_name_html, n_assoc, n_models) {
    list(gene_link, full_name_html, n_assoc, n_models)
  }
)
jsonlite::write_json(list(data = gene_rows), "jekyll/genes.json", auto_unbox = TRUE)

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

panels <- readr::read_tsv("data/panels.tsv", col_types = "cccccc")

panel_counts <- read_tsv(
  "data/all_models.tsv", col_types = cols(PANEL = "c", .default = "-")
) |>
  count(PANEL, name = "N_MODELS")

panels_out <- panels |>
  left_join(panel_counts, by = "PANEL", relationship = "one-to-one")

write_tsv(panels_out, "jekyll/_data/panels.tsv")
