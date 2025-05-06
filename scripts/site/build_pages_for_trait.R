# Load the results of a single trait and generate the basic markdown pages (frontmatter only) for the trait page and each of its hit loci, along with data tables for Jekyll to render in the pages.

# Inputs:
# data/gene_names.tsv
# data/traits.par
# data/panels.par
# data/traits.par.nfo
# data/twas_out/{trait}.all.tsv
# data/twas_out/{trait}.post.tsv

# Outputs:
# jekyll/traits/{trait}.md
# jekyll/traits/{trait}/{locus}.md
# jekyll/_data/trait_loci/{trait}.tsv
# jekyll/_data/trait_pleio/{trait}.tsv
# jekyll/_data/trait_panels/{trait}.tsv

suppressPackageStartupMessages(library(tidyverse))

gene_id <- function(pheno_ID) {
    # Extract gene ID from RNA phenotype ID
    str_extract(pheno_ID, "^[^:.]+")
}

shorten_ids <- function(ID, modality) {
    # Remove extra information from RNA phenotype IDs
    case_when(
        modality %in% c("alternative polyA", "alternative TSS", "isoform ratio") ~ str_extract(ID, "[^:.]+$"),
        modality == "intron excision ratio" ~ str_c("chr", ID |> str_replace("^[^:]+:", "") |> str_replace(":[^:]+$", "")),
        .default = ID
    )
}

arg <- commandArgs(trailingOnly = TRUE)
i_trait <- as.numeric(arg[1])

# Load gene names since Ensembl IDs were used for TWAS
gene_names <- read_tsv("data/gene_names.tsv", col_types = "cc") |>
    deframe()

traits_nfo <- read_tsv("data/traits.par.nfo", col_types = "ciiid")

tbl_traits <- read_tsv("data/traits.par", col_types = "ccicccc") |>
    left_join(traits_nfo, by = "ID")

cat("reading data/panels.par\n")
tbl_panels <- read_tsv("data/panels.par", col_types = "ccccci")

# Count the number of significant genes for the trait
i <- i_trait
trait <- tbl_traits$ID[i]

cat("reading", tbl_traits$OUTPUT[i], "\n")
cur <- read_tsv(
    tbl_traits$OUTPUT[i],
    col_types = cols(PANEL = "c", FILE = "c", ID = "c", TWAS.Z = "d", TWAS.P = "d", .default = "-")
) |>
    filter(!is.na(TWAS.P))

n <- nrow(cur)
top_models <- which(cur$TWAS.P < 0.05 / n)
top_genes <- unique(gene_id(cur$ID[top_models]))

df_cur_models <- cur |>
    summarise(avg.chisq = mean(TWAS.Z^2, na.rm = TRUE),
              num.hits = sum(TWAS.P < 0.05 / n, na.rm = TRUE),
              pct.hits = 100 * mean(TWAS.P < 0.05 / n, na.rm = TRUE),
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

page <- c(
    "---",
    str_glue('title: "{tbl_traits$NAME[i]}"'),
    str_glue('permalink: traits/{trait}/'),
    "layout: trait",
    str_glue("id: {trait}"),
    "---"
)
write_lines(page, fout_page)

# ---- Get clumped and conditional loci
cur_clumps <- read_tsv(
    str_replace(tbl_traits$OUTPUT[i], ".all.tsv", ".post.tsv"), col_types = "ciiiiidddd"
) |>
    arrange(CHR, P0) |>
    mutate(locus_num = seq_len(n()),
           genes = "")

# load clumped genes
clump_mod <- vector()
for (ii in seq_len(nrow(cur_clumps))) {
    cur_genes_tbl <- read_tsv(str_glue("{cur_clumps$FILE[ii]}.genes"),
                              col_types = "ccciiidcdcdddiicdddddddddld") |>
        mutate(num = seq_len(n()),
               gene_id = gene_id(ID),
               link = str_glue("<a href='{{{{ site.baseurl }}}}genes/{gene_id}'>{gene_names[gene_id]}</a>")) |>
    left_join(select(tbl_panels, PANEL, TISSUE, MODALITY), by = "PANEL")
    
    cur_genes <- sort(unique(cur_genes_tbl$gene_id[cur_genes_tbl$JOINT]))
    cur_genes <- str_c(str_glue("<a href='{{{{ site.baseurl }}}}genes/{cur_genes}'>{gene_names[cur_genes]}</a>"), collapse = " ")
    cur_clumps$genes[ii] <- cur_genes
    
    clump_mod <- unique(c(clump_mod, cur_genes_tbl$FILE[cur_genes_tbl$JOINT]))
    
    fout_clump <- str_glue("jekyll/traits/{trait}/{ii}.md")
    cat(str_glue('---\ntitle: "{tbl_traits$NAME[i]}"\npermalink: traits/{trait}/{ii}/\nlayout: locus\n---\n\n'), sep = "", file = fout_clump)
    cat("{: .breadcrumb}\n",
        "[Hub]({{ site.baseurl }}) : [Traits]({{ site.baseurl }}traits) : ",
        tbl_traits$link[i], " : ",
        sep = "", file = fout_clump, append = TRUE)
    if (ii > 1) {
        cat(str_glue(" [ ← ]({{{{ site.baseurl }}}}traits/{trait}/{ii-1}) "), sep = "", file = fout_clump, append = TRUE)
    }
    if (ii < nrow(cur_clumps)) {
        cat(str_glue(" [ → ]({{{{ site.baseurl }}}}traits/{trait}/{ii+1})"), sep = "", file = fout_clump, append = TRUE)
    }
    pos0 <- formatC(cur_clumps$P0[ii], format = "f", big.mark = ",", drop0trailing = TRUE)
    pos1 <- formatC(cur_clumps$P1[ii], format = "f", big.mark = ",", drop0trailing = TRUE)
    cat(str_glue("\n\n# chr{cur_clumps$CHR[ii]}:{pos0}-{pos1}\n\n"), sep = "", file = fout_clump, append = TRUE)
    cat(str_glue("{{}: .text-center}}\nTrait: {tbl_traits$NAME[i]}\n\n"), sep = "", file = fout_clump, append = TRUE)
    cat("{: .text-center}\n",
        "`Best TWAS P=", cur_clumps$BEST.TWAS.P[ii],
        " · Best GWAS P=", cur_clumps$BEST.SNP.P[ii],
        " conditioned to ", cur_clumps$COND.SNP.P[ii],
        "`\n\n", sep = "", file = fout_clump, append = TRUE)
    
    system(str_glue("cp {cur_clumps$FILE[ii]}.cond.csv jekyll/traits/{trait}/{ii}.cond.csv"))
    cat('<div id="graph"></div>\n<script>\nfetch("../', ii, '.cond.csv").then(response => response.text()).then(csvText => { const data = d3.csvParse(csvText); processData(data); });\n</script>\n', sep = "", file = fout_clump, append = TRUE)
    # BCAC.1.post.loc_10.cond.csv
    
    cat("\n### Associated models\n\n", sep = "", file = fout_clump, append = TRUE)
    cat('| # | Tissue | Gene | Modality | RNA phenotype | <span title="Heritability estimate for the given transcriptomic model">h2</span> | eQTL R2 | model | # weights | model R2 | model R2 P | eQTL GWAS Z | TWAS Z | TWAS P | Top SNP corr | <span title="Posterior probability of two distinct causal variants">PP3</span> | <span title="Posterior probability of a single shared causal variant">PP4</span> | <span title="Whether the RNA phenotype is in the joint model">joint</span> |', "| --- |", sep = "\n", file = fout_clump, append = TRUE)
    
    cur_genes_tbl |>
        mutate(COLOC.PP3 = round(COLOC.PP3, 2),
               COLOC.PP4 = round(COLOC.PP4, 2),
               MODELCV.R2 = round(MODELCV.R2, 2),
               HSQ = round(HSQ, 2),
               EQTL.R2 = round(EQTL.R2, 2)) |>
        select(num, TISSUE, link, MODALITY, ID, HSQ, EQTL.R2, MODEL, NWGT, MODELCV.R2, MODELCV.PV, EQTL.GWAS.Z, TWAS.Z, TWAS.P, TOP.SNP.COR, COLOC.PP3, COLOC.PP4, JOINT) |>
        mutate(ID = shorten_ids(ID, MODALITY)) |>
        as.data.frame() |>
        format(digits = 2) |>
        write.table(quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " | ", file = fout_clump, append = TRUE)
    cat(
        "{: #models}\n\n",
        "**h2**: Heritability estimate for the given transcriptomic model. ",
        "**PP3**: Posterior probability of two distinct causal variants. ",
        "**PP4**: Posterior probability of a single shared causal variant. ",
        "**joint**: Whether the RNA phenotype is in the joint model.\n\n",
        file = fout_clump, append = TRUE
    )
}

cur_clumps |>
    mutate(VAR.EXP = round(VAR.EXP * 100, 0)) |>
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
        genes
    ) |>
    mutate(across(where(is.numeric), ~ round(.x, 2))) |>
    write_tsv(fout_loci)

# ---- Get pleiotropic loci
n_pleiot <- 0
for (ii in seq_len(nrow(tbl_traits))) {
    if ((ii != i) && length(top_models) > 0) {
        other_cur <- read_tsv(
            tbl_traits$OUTPUT[ii],
            col_types = cols(FILE = "c", ID = "c", TWAS.Z = "d", TWAS.P = "d", .default = "-")
        ) |>
            filter(!is.na(TWAS.P))
        pleio <- inner_join(cur, other_cur, by = c("FILE", "ID"), relationship = "one-to-one") |>
            filter(TWAS.P.x < 0.05 / n,
                   TWAS.P.y < 0.05 / length(top_models))
        if (nrow(pleio) > 0) {
            if(nrow(pleio) >= 4) {
                tst <- with(pleio, cor.test(TWAS.Z.x, TWAS.Z.y))
            } else {
                tst <- data.frame("est" = 0, "p.value" = 1)
            }
            genes <- sort(unique(gene_id(pleio$ID)))
            genes_link <- str_c(str_glue("<a href='{{{{ site.baseurl }}}}genes/{genes}'>{gene_names[genes]}</a>"), collapse=" ")
            n_genes_twas <- pleio |>
                filter(TWAS.P.y < 0.05 / n) |>
                with(length(unique(gene_id(ID))))
            df_tmp <- data.frame(
                trait_id = tbl_traits$ID[ii],
                chisq_ratio = round(mean(pleio$TWAS.Z.y^2, na.rm = TRUE) / tbl_traits$AVG.CHISQ[ii], 2),
                num_genes = length(genes),
                num_genes_twas = n_genes_twas,
                pct_genes_twas = round(100 * n_genes_twas / tbl_traits$NUM.JOINT.GENES[i], 1), # Replaced 3 with i
                corr = round(tst$est, 2),
                p_val = tst$p.value,
                genes = genes_link
            )
            if (n_pleiot == 0) {
                df_pleiot <- df_tmp
            } else {
                df_pleiot <- rbind(df_pleiot, df_tmp)
            }
            n_pleiot <- nrow(df_pleiot)
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
        genes
    ) |>
    mutate(across(where(is.numeric), ~ round(.x, 2))) |>
    write_tsv(fout_pleio)

df_cur_models |>
    select(
        tissue = TISSUE,
        modality = MODALITY,
        num_hits = num.hits,
        pct_hits = pct.hits,
        avg_chisq = avg.chisq
    ) |>
    mutate(across(where(is.numeric), ~ round(.x, 2))) |>
    write_tsv(fout_panels)
