# Run from the repository root: Rscript tests/test_pleiotropy.R
suppressPackageStartupMessages(library(tidyverse))

test_pleiotropy <- function() {
  root <- getwd()
  renderer <- file.path(root, "scripts/site/build_pleiotropy.R")
  dir.create(".snakemake", showWarnings = FALSE)
  workdir <- tempfile("pleiotropy-test-", tmpdir = file.path(root, ".snakemake"))
  dir.create(file.path(workdir, "data/twas_out"), recursive = TRUE)
  dir.create(file.path(workdir, "jekyll/_data/trait_pleio"), recursive = TRUE)
  setwd(workdir)
  on.exit(setwd(root))

  ids <- c("alpha", "beta", "gamma", "zero", "delta")
  tibble(OUTPUT = str_c("data/twas_out/", ids, ".all.tsv"), ID = ids,
         N = 100, PROJECT = "test", TAGS = "", NAME = ids, DESCRIPTION = "") |>
    write_tsv("data/traits.tsv")
  tibble(ID = ids, NUM.LOCI = c(1, 1, 1, 0, 1),
         NUM.JOINT.GENES = c(3, 1, 3, 0, 1), NUM.GENES = c(3, 2, 3, 0, 1),
         AVG.CHISQ = c(4, 2, 1, 1, 1)) |>
    write_tsv("data/traits.summary.tsv")
  write_lines(c("ID", as.character(1:100)), "data/all_models.tsv")

  # Rows 1 and 6 have the same phenotype ID but different model files.
  keys <- tibble(FILE = str_c("weights/", 1:6),
                 ID = c("GeneA__one", "GeneA__two", "GeneB", "GeneC", "GeneD", "GeneA__one"))
  pvals <- list(c(rep(1e-6, 4), .01, NA), c(1e-6, .005, .0125, 1e-6, .05, NA),
                c(rep(1e-6, 4), .9, NA), rep(.5, 6), c(rep(.5, 5), 1e-6))
  zvals <- list(c(2, 4, 6, 8, 0, NA), c(1, 2, 3, 4, 0, NA),
                c(-1, -2, -3, -4, 0, NA), rep(0, 6), c(rep(0, 5), 2))
  for (i in seq_along(ids)) {
    keys |>
      mutate(TWAS.Z = zvals[[i]], TWAS.P = pvals[[i]]) |>
      write_tsv(str_c("data/twas_out/", ids[i], ".all.tsv"))
  }

  run_renderer <- function() {
    output <- system2("Rscript", shQuote(renderer), stdout = TRUE, stderr = TRUE)
    write_lines(output, "run.log")
    stopifnot(is.null(attr(output, "status")))
  }
  read_result <- function(id) read_tsv(str_c("jekyll/_data/trait_pleio/", id, ".tsv"),
                                      col_types = cols(p_val = "d", .default = col_guess()))
  run_renderer()
  a <- read_result("alpha")
  b <- read_result("beta")
  stopifnot(identical(a$trait_id, c("beta", "gamma")),
            a$num_genes[1] == 2, a$num_genes_twas[1] == 2,
            a$pct_genes_twas[1] == 66.7, a$chisq_ratio[1] == 3.5,
            a$corr[1] == 0, a$p_val[1] == 1, a$genes[1] == "GeneA,GeneC",
            a$num_genes[2] == 3, a$corr[2] == -1, a$p_val[2] < 1e-6,
            b$chisq_ratio[b$trait_id == "alpha"] == 8.5,
            # Preserve the original denominator (joint genes), even above 100%.
            b$pct_genes_twas[b$trait_id == "alpha"] == 200,
            nrow(read_result("zero")) == 0, nrow(read_result("delta")) == 0)

  # An empty union is valid when every trait genuinely has no significant models.
  for (id in ids) {
    keys |> mutate(TWAS.Z = 0, TWAS.P = .5) |>
      write_tsv(str_c("data/twas_out/", id, ".all.tsv"))
  }
  run_renderer()
  stopifnot(all(vapply(ids, function(id) nrow(read_result(id)) == 0, logical(1))))
  cat("Pleiotropy statistical and empty-union checks passed. Evidence:", workdir, "\n")
}

test_pleiotropy()
