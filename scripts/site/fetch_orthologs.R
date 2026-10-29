# Fetch the model-gene ortholog mapping from Ensembl release 108.
# Association counts are added separately by build_cross_species.R.
library(biomaRt)
suppressPackageStartupMessages(library(tidyverse))

human_genes <- read_delim("data/cross_species/genes.models.nfo",
                        col_names = c("human_symbol", "count"), col_types = "ci")
rat_models <- read_tsv("data/all_models.tsv",
                       col_types = cols(ID = "c", .default = "-")) |>
  transmute(rat_symbol = str_replace(ID, "__.+$", "")) |>
  distinct()

# Set up biomaRt connections
mart_human <- useEnsembl(
  biomart = "ensembl",
  dataset = "hsapiens_gene_ensembl",
  host = "https://oct2022.archive.ensembl.org", # v108
  mirror = "useast"
)
mart_rat <- useEnsembl(
  biomart = "ensembl",
  dataset = "rnorvegicus_gene_ensembl",
  host = "https://oct2022.archive.ensembl.org", # v108
  mirror = "useast"
)

# Get human Ensembl IDs from HGNC symbols
human_names <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "hgnc_symbol",
  values = human_genes$human_symbol,
  mart = mart_human
) |> 
  rename(
    human_id = ensembl_gene_id,
    human_symbol = hgnc_symbol
  )

# Get rat Ensembl IDs from symbols
rat_names <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "external_gene_name",
  values = unique(rat_models$rat_symbol),
  mart = mart_rat
) |>
  as_tibble() |>
  rename(rat_id = ensembl_gene_id, rat_symbol = external_gene_name) |>
  distinct()

# Map orthologs in both directions so all genes in both lists are included
human_to_rat <- getBM(
  attributes = c("ensembl_gene_id", "rnorvegicus_homolog_ensembl_gene", "rnorvegicus_homolog_associated_gene_name"),
  filters = "ensembl_gene_id",
  values = human_names$human_id,
  mart = mart_human
) |>
  as_tibble() |>
  rename(
    human_id = ensembl_gene_id,
    rat_id = rnorvegicus_homolog_ensembl_gene,
    rat_symbol = rnorvegicus_homolog_associated_gene_name,
  ) |>
  # There are a few IDs listed more than once
  left_join(human_names, by = "human_id", relationship = "many-to-many") |>
  select(human_id, human_symbol, rat_id, rat_symbol)

rat_to_human <- getBM(
  attributes = c("ensembl_gene_id", "hsapiens_homolog_ensembl_gene", "hsapiens_homolog_associated_gene_name"),
  filters = "ensembl_gene_id",
  values = unique(rat_names$rat_id[!is.na(rat_names$rat_id)]),
  mart = mart_rat
) |>
  as_tibble() |>
  rename(
    rat_id = ensembl_gene_id,
    human_id = hsapiens_homolog_ensembl_gene,
    human_symbol = hsapiens_homolog_associated_gene_name,
  ) |>
  left_join(rat_names, by = "rat_id", relationship = "many-to-one") |>
  select(human_id, human_symbol, rat_id, rat_symbol)

orthologs <- bind_rows(human_to_rat, rat_to_human) |>
  distinct() |>
  mutate(
    human_id = if_else(human_id == "", NA, human_id),
    human_symbol = if_else(human_symbol == "", NA, human_symbol),
    rat_id = if_else(rat_id == "", NA, rat_id),
    rat_symbol = if_else(rat_symbol == "", NA, rat_symbol),
  ) |>
  arrange(human_symbol, rat_symbol)

write_tsv(orthologs, "data/cross_species/human_rat_orthologs.ensembl108.tsv")
