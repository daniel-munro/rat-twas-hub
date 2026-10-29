# Render one locus from the trait's chromosome/position-sorted post report.
suppressPackageStartupMessages(library(tidyverse))

pheno_to_gene_id <- function(pheno_id) str_replace(pheno_id, "__.+$", "")

shorten_ids <- function(pheno_id, modality) {
  # Remove extra information from RNA phenotype IDs
  # alternative TSS/polyA:            {gene_id}__grp_{grp#}_{up|down}stream_{transcript_id}      -> {transcript_id}
  # gene expression & mRNA stability: {gene_id}                                                  -> {gene_id}
  # intron excision ratio:            {gene_id}__{chrom}_{pos1}_{pos2}_clu_{cluster_id}_{strand} -> {chrom}_{pos1}_{pos2}
  # isoform ratio:                    {gene_id}__{transcript_id}                                 -> {transcript_id}
  dplyr::case_when(
    modality %in% c("alternative TSS", "alternative polyA") ~ pheno_id |> stringr::str_replace("^.+stream_", ""),
    modality == "intron excision ratio" ~ pheno_id |> stringr::str_replace("^.+__", "") |> stringr::str_replace("_clu_.+$", ""),
    modality == "isoform ratio" ~ pheno_id |> stringr::str_replace("^.+__", ""),
    .default = pheno_id
  )
}

args <- commandArgs(trailingOnly = TRUE)
trait <- args[[1]]
i_locus <- as.integer(args[[2]])
tbl_traits <- read_tsv("data/traits.tsv", col_types = "ccicccc", num_threads = 1)
i_trait <- match(trait, tbl_traits$ID)
tbl_panels <- read_tsv("data/panels.tsv", col_types = "ccccci", num_threads = 1)
trait_loci <- read_tsv(str_glue("data/twas_out/{trait}.post.tsv"), col_types = "cciiiidddd", num_threads = 1) |>
  mutate(chr_num = as.integer(str_replace(CHR, "chr", ""))) |>
  arrange(chr_num, P0)
stopifnot(!is.na(i_trait), i_locus >= 1L, i_locus <= nrow(trait_loci))

print(i_locus)
fout_locus_page <- str_glue("jekyll/traits/{trait}/{i_locus}.md")
fout_locus_models <- str_glue("jekyll/_data/locus_models/{trait}/{i_locus}.tsv")

stopifnot(file.copy(
  str_glue("{trait_loci$FILE[i_locus]}.cond.csv"),
  str_glue("jekyll/traits/{trait}/{i_locus}.cond.csv"), overwrite = TRUE
))

locus_models <- read_tsv(
  str_glue("{trait_loci$FILE[i_locus]}.genes"),
  col_types = "cccciidcdcdddiicdddddddddld", num_threads = 1
) |>
  mutate(num = seq_len(n()),
         gene_id = pheno_to_gene_id(ID)) |>
  left_join(select(tbl_panels, PANEL, TISSUE, MODALITY), by = "PANEL")

pos0 <- formatC(trait_loci$P0[i_locus], format = "f", big.mark = ",", drop0trailing = TRUE)
pos1 <- formatC(trait_loci$P1[i_locus], format = "f", big.mark = ",", drop0trailing = TRUE)
c(
  "---",
  str_glue('title: "{tbl_traits$NAME[i_trait]}"'),
  str_glue("permalink: traits/{trait}/{i_locus}/"),
  "layout: locus",
  str_glue("trait_id: {trait}"),
  str_glue("locus_num: {i_locus}"),
  str_glue('pos0: "{pos0}"'),
  str_glue('pos1: "{pos1}"'),
  "---"
) |>
  write_lines(fout_locus_page)

locus_models |>
  select(
    num,
    tissue = TISSUE,
    gene_id,
    modality = MODALITY,
    id = ID,
    hsq = HSQ,
    eqtl_r2 = EQTL.R2,
    model = MODEL,
    nwgt = NWGT,
    modelcv_r2 = MODELCV.R2,
    modelcv_pv = MODELCV.PV,
    eqtl_gwas_z = EQTL.GWAS.Z,
    twas_z = TWAS.Z,
    twas_p = TWAS.P,
    top_snp_corr = TOP.SNP.COR,
    coloc_pp3 = COLOC.PP3,
    coloc_pp4 = COLOC.PP4,
    joint = JOINT,
  ) |>
  mutate(
    id = shorten_ids(id, modality),
    hsq = round(hsq, 2),
    eqtl_r2 = round(eqtl_r2, 2),
    modelcv_r2 = round(modelcv_r2, 2),
    modelcv_pv = format(modelcv_pv, digits = 3, scientific = TRUE),
    eqtl_gwas_z = round(eqtl_gwas_z, 2),
    twas_z = round(twas_z, 2),
    twas_p = format(twas_p, digits = 3, scientific = TRUE),
    top_snp_corr = round(top_snp_corr, 2),
    coloc_pp3 = round(coloc_pp3, 2),
    coloc_pp4 = round(coloc_pp4, 2),
  ) |>
  write_tsv(fout_locus_models)
