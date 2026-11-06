# Generate a trait page and summary tables; locus pages are rendered separately.

# Inputs:
# data/traits.tsv
# data/panels.tsv
# data/twas_out/{trait}.all.tsv
# data/twas_out/{trait}.post.tsv
# data/twas_out/{trait}/{trait}.{chrom}.post.loc_{locus}.genes

# Outputs:
# jekyll/traits/{trait}.md
# jekyll/_data/trait_loci/{trait}.tsv
# jekyll/_data/trait_panels/{trait}.tsv
# jekyll/_data/trait_hits/{trait}.tsv

suppressPackageStartupMessages(library(tidyverse))

pheno_to_gene_id <- function(pheno_id) {
  # Extract gene ID from RNA phenotype ID
  stringr::str_replace(pheno_id, "__.+$", "")
}

arg <- commandArgs(trailingOnly = TRUE)
trait <- as.character(arg[1])

tbl_traits <- read_tsv("data/traits.tsv", col_types = "ccicccc")

# Determine the index of the requested trait_id and validate
i_trait <- match(trait, tbl_traits$ID)
if (is.na(i_trait)) {
  stop(stringr::str_glue("Trait ID '{trait}' not found in data/traits.tsv"))
}

tbl_panels <- read_tsv("data/panels.tsv", col_types = "ccccci")

cat("reading", tbl_traits$OUTPUT[i_trait], "\n")
# Preserve the source representation of very small P-values in the hit table.
trait_assocs <- read_tsv(
  tbl_traits$OUTPUT[i_trait],
  col_types = cols(PANEL = "c", FILE = "c", ID = "c", TWAS.Z = "d", TWAS.P = "c", .default = "-")
) |>
  mutate(p = as.numeric(TWAS.P)) |>
  filter(!is.na(p)) |>
  mutate(gene_id = pheno_to_gene_id(ID))

# Match TWAS hit selection and plots, including models without a result.
n_models <- length(read_lines("data/all_models.tsv")) - 1
threshold <- 0.05 / n_models

trait_panels <- trait_assocs |>
  summarise(hit_genes = n_distinct(gene_id[p < threshold]),
            num_hits = sum(p < threshold),
            num_tested = n(),
            .by = PANEL) |>
  left_join(tbl_panels, by = "PANEL")

# ---- PRINT TRAIT PAGE

fout_page <- str_glue("jekyll/traits/{trait}.md")
fout_loci <- str_glue("jekyll/_data/trait_loci/{trait}.tsv")
fout_panels <- str_glue("jekyll/_data/trait_panels/{trait}.tsv")

# These directories contain only this trait's generated locus files. Clearing
# them removes obsolete loci before the downstream locus jobs recreate outputs.
unlink(str_glue("jekyll/traits/{trait}"), recursive = TRUE)
unlink(str_glue("jekyll/_data/locus_models/{trait}"), recursive = TRUE)

c(
  "---",
  str_glue('title: "{tbl_traits$NAME[i_trait]}"'),
  str_glue("permalink: traits/{trait}/"),
  "layout: trait",
  str_glue("id: {trait}"),
  "---"
) |>
  write_lines(fout_page)

# ---- Get clumped and conditional loci
trait_loci <- read_tsv(
  str_replace(tbl_traits$OUTPUT[i_trait], ".all.tsv", ".post.tsv"), col_types = "cciiiicccd"
) |>
  mutate(chr_num = as.integer(str_replace(CHR, "chr", ""))) |>
  arrange(chr_num, P0) |>
  select(-chr_num) |>
  mutate(locus_num = seq_len(n()),
         genes = "")

# Joint genes are needed for the trait's locus table; locus pages are separate
# tracked outputs rendered by build_locus_page.R.
locus_keys <- vector("list", nrow(trait_loci))
for (i_locus in seq_len(nrow(trait_loci))) {
  locus_models <- read_tsv(
    str_glue("{trait_loci$FILE[i_locus]}.genes"),
    col_types = cols(PANEL = "c", FILE = "c", ID = "c", JOINT = "l", .default = "-")
  )
  trait_loci$genes[i_locus] <- str_c(
    sort(unique(pheno_to_gene_id(locus_models$ID[locus_models$JOINT]))), collapse = ","
  )
  locus_keys[[i_locus]] <- distinct(locus_models, PANEL, FILE, ID) |>
    mutate(locus = i_locus)
}

# Match exact model rows rather than assigning a hit by the gene's coordinates.
# Some significant models have no post-processing locus.
hit_loci <- bind_rows(tibble(PANEL = character(), FILE = character(), ID = character(),
                             locus = integer()), bind_rows(locus_keys)) |>
  summarise(locus = paste(sort(unique(locus)), collapse = ","), .by = c(PANEL, FILE, ID))
trait_hits <- trait_assocs |>
  filter(p < threshold) |>
  left_join(hit_loci, by = c("PANEL", "FILE", "ID"), relationship = "many-to-one") |>
  left_join(select(tbl_panels, PANEL, TISSUE, MODALITY),
            by = "PANEL", relationship = "many-to-one") |>
  arrange(p, desc(abs(TWAS.Z)), PANEL, ID, FILE) |>
  select(gene_id, tissue = TISSUE, modality = MODALITY, pheno_id = ID,
         twas_z = TWAS.Z, twas_p = TWAS.P, locus)
dir.create("jekyll/_data/trait_hits", showWarnings = FALSE, recursive = TRUE)
write_tsv(trait_hits, str_glue("jekyll/_data/trait_hits/{trait}.tsv"), na = "")

trait_loci |>
  select(
    locus_num,
    chr = CHR,
    p0 = P0,
    p1 = P1,
    hit_genes = HIT.GENES,
    joint_genes = JOINT.GENES,
    best_twas_p = BEST.TWAS.P,
    best_snp_p = BEST.SNP.P,
    cond_snp_p = COND.SNP.P,
    var_exp = VAR.EXP,
    genes,
  ) |>
  mutate(
    across(c(best_twas_p, best_snp_p, cond_snp_p),
           ~ format(as.numeric(.x), digits = 3, scientific = TRUE)),
    var_exp = round(var_exp * 100, 0),
  ) |>
  write_tsv(fout_loci)

trait_panels |>
  select(
    tissue = TISSUE,
    modality = MODALITY,
    hit_genes,
    num_hits,
    num_tested,
  ) |>
  write_tsv(fout_panels)
