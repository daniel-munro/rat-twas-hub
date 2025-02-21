suppressPackageStartupMessages(library(tidyverse))

sample_size <- read_tsv("data/gwas_original/sample_sizes.tsv",
                        col_names = c("project_id", "trait_id", "N"),
                        col_types = "ccc")

df <- read_tsv("data/gwas_original/traits.tsv", col_types = "ccccc") |>
  left_join(sample_size, by = c("project_id", "trait_id"), relationship = "one-to-one") |>
  mutate(OUTPUT = str_glue("data/twas_out/{trait_id}.all.tsv")) |>
  select(OUTPUT,
         ID = trait_id,
         N,
         PROJECT = project_id,
         TAGS = tags,
         NAME = name,
         DESCRIPTION = description)

write_tsv(df, "data/traits.par")
