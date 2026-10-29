# Build cross-species ortholog table

# Inputs:
# data/cross_species/genes.models.nfo    human: gene symbol and model count (space-separated)
# data/genes_n_models.tsv                rat: gene symbol and model count (tab-separated)
# data/cross_species/genes.nfo           human: gene symbol and association count (space-separated)
# data/genes_n_assoc.tsv                 rat: gene symbol and association count (tab-separated)

# data/cross_species/human_rat_orthologs.ensembl108.tsv  frozen mapping

# Outputs:
# jekyll/cross-species.json

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

# Reuse the frozen mapping; association updates do not query Ensembl.
orthologs <- read_tsv(
  "data/cross_species/human_rat_orthologs.ensembl108.tsv", col_types = "cccc"
)

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
