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

  ids <- c("alpha", "beta", "gamma", "zero", "delta", "epsilon", "eta")
  tibble(OUTPUT = str_c("data/twas_out/", ids, ".all.tsv"), ID = ids,
         N = 100, PROJECT = "test", TAGS = "", NAME = ids, DESCRIPTION = "") |>
    write_tsv("data/traits.tsv")
  tibble(ID = ids, NUM.LOCI = 1, NUM.JOINT.GENES = 1, NUM.GENES = 3,
         AVG.CHISQ = c(4, 2, 1, 1, 1, 1, 1)) |>
    write_tsv("data/traits.summary.tsv")
  write_lines(c("ID", as.character(1:100)), "data/all_models.tsv")

  # Rows 1 and 6 have the same phenotype ID but different model files.
  keys <- tibble(FILE = str_c("weights/", 1:7),
                 ID = c("GeneA__one", "GeneA__two", "GeneB", "GeneC", "GeneD", "GeneA__one", "GeneE"))
  pvals <- list(
    alpha = c(rep(1e-6, 4), .01, NA, 1e-6),
    beta = c(1e-6, .005, .01, 1e-6, .05, NA, .005),
    gamma = c(rep(1e-6, 4), .9, NA, .0005),
    zero = rep(.5, 7),
    delta = c(rep(.5, 5), 1e-6, .5),
    epsilon = c(rep(.005, 4), rep(.5, 3)),
    eta = c(rep(1e-6, 4), rep(.5, 3))
  )
  zvals <- list(
    alpha = c(2, 4, 6, 8, 0, NA, 9),
    beta = c(1, 2, 3, 4, 0, NA, NA),
    gamma = c(-1, -2, -3, -4, 0, NA, -4.5),
    zero = rep(0, 7),
    delta = c(rep(0, 5), 2, 0),
    epsilon = rep(2, 7),
    eta = c(1, -1, -1, 1, 0, 0, 0)
  )
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
  alpha <- read_result("alpha")
  stopifnot(identical(alpha$trait_id, c("beta", "gamma", "epsilon", "eta")),
            !("pct_genes_twas" %in% names(alpha)))

  beta <- filter(alpha, trait_id == "beta")
  stopifnot(beta$num_genes_twas == 2, beta$num_models_twas == 2,
            beta$genes_twas == "GeneA,GeneC",
            beta$num_genes == 3, beta$genes == "GeneA,GeneC,GeneE",
            beta$chisq_ratio == 3.5, is.na(beta$corr), is.na(beta$p_val))
  reverse <- filter(read_result("beta"), trait_id == "alpha")
  stopifnot(reverse$num_models_twas == beta$num_models_twas,
            reverse$genes_twas == beta$genes_twas, reverse$num_genes == 2,
            reverse$chisq_ratio == 8.5)

  gamma <- filter(alpha, trait_id == "gamma")
  stopifnot(gamma$num_genes_twas == 3, gamma$num_models_twas == 4,
            gamma$genes_twas == "GeneA,GeneB,GeneC",
            gamma$corr == -1, gamma$p_val < 1e-6)
  epsilon <- filter(alpha, trait_id == "epsilon")
  stopifnot(epsilon$num_genes_twas == 0, epsilon$num_models_twas == 0,
            epsilon$num_genes == 3, is.na(epsilon$corr), is.na(epsilon$p_val))
  eta <- filter(alpha, trait_id == "eta")
  stopifnot(eta$corr == 0, eta$p_val == 1)
  stopifnot(nrow(read_result("zero")) == 0, nrow(read_result("delta")) == 0,
            nrow(read_result("epsilon")) == 0)

  # An empty union is valid when every trait genuinely has no significant models.
  for (id in ids) {
    keys |> mutate(TWAS.Z = 0, TWAS.P = .5) |>
      write_tsv(str_c("data/twas_out/", id, ".all.tsv"))
  }
  run_renderer()
  stopifnot(all(vapply(ids, function(id) nrow(read_result(id)) == 0, logical(1))))
  cat("Shared xTWAS and directional test checks passed. Evidence:", workdir, "\n")
}

test_pleiotropy()
