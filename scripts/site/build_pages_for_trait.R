# Generate a trait page and summary tables; locus pages are rendered separately.

# Inputs:
# data/traits.tsv
# data/panels.tsv
# data/traits.summary.tsv
# data/twas_out/{trait}.all.tsv
# data/twas_out/{trait}.post.tsv
# data/twas_out/{trait}/{trait}.{chrom}.post.loc_{locus}.cond.csv

# Outputs:
# jekyll/traits/{trait}.md
# jekyll/_data/trait_loci/{trait}.tsv
# jekyll/_data/trait_pleio/{trait}.tsv
# jekyll/_data/trait_panels/{trait}.tsv

suppressPackageStartupMessages(library(tidyverse))

pheno_to_gene_id <- function(pheno_id) {
  # Extract gene ID from RNA phenotype ID
  stringr::str_replace(pheno_id, "__.+$", "")
}

arg <- commandArgs(trailingOnly = TRUE)
trait <- as.character(arg[1])

traits_nfo <- read_tsv("data/traits.summary.tsv", col_types = "ciiid")

tbl_traits <- read_tsv("data/traits.tsv", col_types = "ccicccc") |>
  left_join(traits_nfo, by = "ID")

# Determine the index of the requested trait_id and validate
i_trait <- match(trait, tbl_traits$ID)
if (is.na(i_trait)) {
  stop(stringr::str_glue("Trait ID '{trait}' not found in data/traits.tsv"))
}

tbl_panels <- read_tsv("data/panels.tsv", col_types = "ccccci")

cat("reading", tbl_traits$OUTPUT[i_trait], "\n")
trait_assocs <- read_tsv(
  tbl_traits$OUTPUT[i_trait],
  col_types = cols(PANEL = "c", FILE = "c", ID = "c", TWAS.Z = "d", TWAS.P = "d", .default = "-")
) |>
  filter(!is.na(TWAS.P))

# Match TWAS hit selection and plots, including models without a result.
n_models <- length(read_lines("data/all_models.tsv")) - 1
top_models <- which(trait_assocs$TWAS.P < 0.05 / n_models)
top_genes <- unique(pheno_to_gene_id(trait_assocs$ID[top_models]))

trait_panels <- trait_assocs |>
  summarise(avg_chisq = mean(TWAS.Z^2, na.rm = TRUE),
            num_hits = sum(TWAS.P < 0.05 / n_models, na.rm = TRUE),
            pct_hits = 100 * mean(TWAS.P < 0.05 / n_models, na.rm = TRUE),
            .by = PANEL) |>
  left_join(tbl_panels, by = "PANEL")

# ---- PRINT TRAIT PAGE

fout_page <- str_glue("jekyll/traits/{trait}.md")
fout_loci <- str_glue("jekyll/_data/trait_loci/{trait}.tsv")
fout_pleio <- str_glue("jekyll/_data/trait_pleio/{trait}.tsv")
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
  str_replace(tbl_traits$OUTPUT[i_trait], ".all.tsv", ".post.tsv"), col_types = "cciiiidddd"
) |>
  mutate(chr_num = as.integer(str_replace(CHR, "chr", ""))) |>
  arrange(chr_num, P0) |>
  select(-chr_num) |>
  mutate(locus_num = seq_len(n()),
         genes = "")

# Joint genes are needed for the trait's locus table; locus pages are separate
# tracked outputs rendered by build_locus_page.R.
for (i_locus in seq_len(nrow(trait_loci))) {
  locus_models <- read_tsv(
    str_glue("{trait_loci$FILE[i_locus]}.genes"),
    col_types = cols(ID = "c", JOINT = "l", .default = "-")
  )
  trait_loci$genes[i_locus] <- str_c(
    sort(unique(pheno_to_gene_id(locus_models$ID[locus_models$JOINT]))), collapse = ","
  )
}

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
    best_twas_p = format(best_twas_p, digits = 3, scientific = TRUE),
    best_snp_p = format(best_snp_p, digits = 3, scientific = TRUE),
    cond_snp_p = format(cond_snp_p, digits = 3, scientific = TRUE),
    var_exp = round(var_exp * 100, 0),
  ) |>
  write_tsv(fout_loci)

# ---- Get pleiotropic associations
df_pleiot <- tibble(
  trait_id = character(),
  chisq_ratio = numeric(),
  num_genes = integer(),
  num_genes_twas = integer(),
  pct_genes_twas = numeric(),
  corr = numeric(),
  p_val = numeric(),
  genes = character()
)

for (i_trait2 in seq_len(nrow(tbl_traits))) {
  print(i_trait2)
  if ((i_trait2 != i_trait) && length(top_models) > 0) {
    trait2_assocs <- read_tsv(
      tbl_traits$OUTPUT[i_trait2],
      col_types = cols(FILE = "c", ID = "c", TWAS.Z = "d", TWAS.P = "d", .default = "-")
    ) |>
      filter(!is.na(TWAS.P))
    pleio <- inner_join(trait_assocs, trait2_assocs, by = c("FILE", "ID"), relationship = "one-to-one") |>
      filter(TWAS.P.x < 0.05 / n_models,
             TWAS.P.y < 0.05 / length(top_models))
    if (nrow(pleio) > 0) {
      if(nrow(pleio) >= 4) {
        tst <- with(pleio, cor.test(TWAS.Z.x, TWAS.Z.y))
      } else {
        tst <- tibble(est = 0, p.value = 1)
      }
      genes <- sort(unique(pheno_to_gene_id(pleio$ID)))
      n_genes_twas <- pleio |>
        filter(TWAS.P.y < 0.05 / n_models) |>
        with(length(unique(pheno_to_gene_id(ID))))
      df_pleiot <- bind_rows(df_pleiot, tibble(
        trait_id = tbl_traits$ID[i_trait2],
        chisq_ratio = mean(pleio$TWAS.Z.y^2, na.rm = TRUE) / tbl_traits$AVG.CHISQ[i_trait2],
        num_genes = length(genes),
        num_genes_twas = n_genes_twas,
        pct_genes_twas = round(100 * n_genes_twas / tbl_traits$NUM.JOINT.GENES[i_trait], 1),
        corr = tst$est,
        p_val = tst$p.value,
        genes = str_c(genes, collapse = ",")
      ))
    }
  }
}

df_pleiot |>
  replace_na(list(pct_genes_twas = 0)) |>
  select(
    trait_id,
    chisq_ratio,
    num_genes,
    num_genes_twas,
    pct_genes_twas,
    corr,
    p_val,
    genes,
  ) |>
  mutate(
    chisq_ratio = round(chisq_ratio, 2),
    corr = round(corr, 2),
    p_val = format(p_val, digits = 3, scientific = TRUE),
  ) |>
  write_tsv(fout_pleio)

trait_panels |>
  select(
    tissue = TISSUE,
    modality = MODALITY,
    num_hits,
    pct_hits,
    avg_chisq,
  ) |>
  mutate(
    pct_hits = round(pct_hits, 1),
    avg_chisq = round(avg_chisq, 2),
  ) |>
  write_tsv(fout_panels)
