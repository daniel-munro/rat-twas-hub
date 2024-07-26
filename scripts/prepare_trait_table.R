suppressPackageStartupMessages(library(tidyverse))

descs <- read_tsv("data/GWAS_original/trait_info.tsv", col_types = "ccc")

sample_size <- read_tsv("data/GWAS_original/sample_sizes.tsv",
                        col_names = c("project", "trait", "N"),
                        col_types = "ccc")

df <- read_tsv("data/GWAS_original/traits_pruned.r2_50.tsv",
               col_names = c("project", "trait"),
               col_types = "cc") |>
  left_join(sample_size, by = c("project", "trait"), relationship = "one-to-one") |>
  left_join(descs, by = c("project", "trait"), relationship = "one-to-one") |>
  mutate(OUTPUT = str_glue("data/twas_out/{trait}.dat"),
         YEAR = "2024",
         TYPE = project) |>
  select(OUTPUT, ID = trait, N, REF = project, YEAR, NAME = trait_description, TYPE)

write_tsv(df, "data/traits.par")
