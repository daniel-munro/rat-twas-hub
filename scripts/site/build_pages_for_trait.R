# Load the results of a single trait and generate the basic markdown pages (frontmatter only) for the trait page and each of its hit loci, along with data tables for Jekyll to render in the pages.

# Inputs:
# data/traits.par
# data/panels.par
# data/traits.par.nfo
# data/twas_out/{trait}.all.tsv
# data/twas_out/{trait}.post.tsv
# data/twas_out/{trait}/{trait}.{chrom}.post.loc_{locus}.cond.csv

# Outputs:
# jekyll/traits/{trait}.md
# jekyll/traits/{trait}/{locus}.md
# jekyll/traits/{trait}/{locus}.cond.csv
# jekyll/_data/trait_loci/{trait}.tsv
# jekyll/_data/trait_pleio/{trait}.tsv
# jekyll/_data/trait_panels/{trait}.tsv
# jekyll/_data/locus_models/{trait}/{locus}.tsv

suppressPackageStartupMessages(library(tidyverse))

pheno_to_gene_id <- function(pheno_id) {
    # Extract gene ID from RNA phenotype ID
    str_extract(pheno_id, "^[^:.]+")
}

shorten_ids <- function(id, modality) {
    # Remove extra information from RNA phenotype IDs
    case_when(
        modality %in% c("alternative polyA", "alternative TSS", "isoform ratio") ~ str_extract(id, "[^:.]+$"),
        modality == "intron excision ratio" ~ str_c("chr", id |> str_replace("^[^:]+:", "") |> str_replace(":[^:]+$", "")),
        .default = id
    )
}

arg <- commandArgs(trailingOnly = TRUE)
i_trait <- as.numeric(arg[1])

traits_nfo <- read_tsv("data/traits.par.nfo", col_types = "ciiid")

tbl_traits <- read_tsv("data/traits.par", col_types = "ccicccc") |>
    left_join(traits_nfo, by = "ID")

tbl_panels <- read_tsv("data/panels.par", col_types = "ccccci")

# Count the number of significant genes for the trait
trait <- tbl_traits$ID[i_trait]

cat("reading", tbl_traits$OUTPUT[i_trait], "\n")
trait_assocs <- read_tsv(
    tbl_traits$OUTPUT[i_trait],
    col_types = cols(PANEL = "c", FILE = "c", ID = "c", TWAS.Z = "d", TWAS.P = "d", .default = "-")
) |>
    filter(!is.na(TWAS.P))

n_models <- nrow(trait_assocs)
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

system(str_glue("mkdir -p jekyll/traits/{trait}"))
system(str_glue("mkdir -p jekyll/_data/trait_loci"))
system(str_glue("mkdir -p jekyll/_data/trait_pleio"))
system(str_glue("mkdir -p jekyll/_data/trait_panels"))

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
    str_replace(tbl_traits$OUTPUT[i_trait], ".all.tsv", ".post.tsv"), col_types = "ciiiiidddd"
) |>
    arrange(CHR, P0) |>
    mutate(locus_num = seq_len(n()),
           genes = "")

for (i_locus in seq_len(nrow(trait_loci))) {
    print(i_locus)
    fout_locus_page <- str_glue("jekyll/traits/{trait}/{i_locus}.md")
    fout_locus_models <- str_glue("jekyll/_data/locus_models/{trait}/{i_locus}.tsv")

    system(str_glue("mkdir -p jekyll/_data/locus_models/{trait}"))
    system(str_glue("cp {trait_loci$FILE[i_locus]}.cond.csv jekyll/traits/{trait}/{i_locus}.cond.csv"))

    locus_models <- read_tsv(
        str_glue("{trait_loci$FILE[i_locus]}.genes"),
        col_types = "ccciiidcdcdddiicdddddddddld"
    ) |>
        mutate(num = seq_len(n()),
               gene_id = pheno_to_gene_id(ID)) |>
    left_join(select(tbl_panels, PANEL, TISSUE, MODALITY), by = "PANEL")
    
    locus_joint_genes <- sort(unique(locus_models$gene_id[locus_models$JOINT]))
    trait_loci$genes[i_locus] <- str_c(locus_joint_genes, collapse = ",")
    
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
    # print(i_trait2)
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
