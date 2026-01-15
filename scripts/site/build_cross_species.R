# Build cross-species ortholog table

# Inputs:
# data/cross_species/genes.models.nfo    human: gene symbol and model count (space-separated)
# data/genes_n_models.tsv                rat: gene symbol and model count (tab-separated)
# data/cross_species/genes.nfo           human: gene symbol and association count (space-separated)
# data/genes_n_assoc.tsv                 rat: gene symbol and association count (tab-separated)

# Outputs:
# jekyll/cross-species.json

library(biomaRt)
suppressPackageStartupMessages(library(tidyverse))

# Read input files
human_genes <- read_delim(
  "data/cross_species/genes.models.nfo",
  col_names = c("human_symbol", "count"),
  col_types = "ci"
)

# Rat inputs are now gene symbols end-to-end
rat_models <- read_delim(
  "data/genes_n_models.tsv",
  col_names = c("rat_symbol", "count"),
  col_types = "ci",
  delim = "\t"
)

human_assoc <- read_delim(
  "data/cross_species/genes.nfo",
  col_names = c("human_symbol", "human_assoc"),
  col_types = "ci"
)

rat_assoc <- read_delim(
  "data/genes_n_assoc.tsv",
  col_names = c("rat_symbol", "rat_assoc"),
  col_types = "ci",
  delim = "\t"
)

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

# Record which human and rat genes list have models
orthologs <- orthologs |>
  mutate(
    human_has_model = human_symbol %in% human_genes$human_symbol,
    rat_has_model = rat_symbol %in% rat_models$rat_symbol
  )

# Record how many associations each human and rat gene has
orthologs <- orthologs |>
  left_join(human_assoc, by = "human_symbol", relationship = "many-to-one") |>
  mutate(human_assoc = if_else(is.na(human_assoc), 0, human_assoc)) |>
  left_join(rat_assoc, by = "rat_symbol", relationship = "many-to-one") |>
  mutate(rat_assoc = if_else(is.na(rat_assoc), 0, rat_assoc))

# write_tsv(orthologs, "data/cross_species/human_rat_orthologs.tsv")
# write_json(
#     list(data = orthologs |> as.list()),
#     "data/cross_species/human_rat_orthologs.json"
# )

orthologs <- orthologs |>
  replace_na(list(human_symbol = "", human_id = "", rat_symbol = "", rat_id = "")) |>
  mutate(
    human_symbol_link = if_else(
      human_has_model,
      str_glue('<em><a href=\\"http://twas-hub.org/genes/{human_symbol}\\">{human_symbol}</a></em>'),
      human_symbol
    ),
    human_id_link = if_else(
      human_has_model,
      str_glue('<em><a href=\\"http://twas-hub.org/genes/{human_symbol}\\">{human_id}</a></em>'),
      human_id
    ),
    rat_symbol_link = if_else(
      rat_has_model,
      str_glue('<em><a href=\\"../genes/{rat_symbol}\\">{rat_symbol}</a></em>'),
      rat_symbol
    ),
    rat_id_link = if_else(
      rat_has_model,
      str_glue('<em><a href=\\"../genes/{rat_symbol}\\">{rat_id}</a></em>'),
      rat_id
    ),
    line = str_glue('["{human_symbol_link}","{human_id_link}","{rat_symbol_link}","{rat_id_link}",{human_assoc},{rat_assoc}]')
  )

cat(
  '{\n"data":[\n',
  str_c(orthologs$line, collapse = ",\n"),
  "]\n}",
  sep = "",
  file = "jekyll/cross-species.json"
)
