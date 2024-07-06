# Load the results of a single trait and generate a markdown file for the trait page and each of its hit loci

# Inputs:
# data/gene_names.tsv
# data/traits.par
# data/panels.par
# data/traits.par.nfo
# data/twas_out/{trait}.dat
# data/twas_out/{trait}.dat.post.report

# Outputs:
# jekyll/traits/{trait}.md
# jekyll/traits/{trait}/{locus}.md

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

tbl_traits <- read_tsv("data/traits.par", col_types = "ccicicc") |>
    mutate(link = str_glue("[{NAME}]({{{{ site.baseurl }}}}traits/{ID})")) |>
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

fout <- str_glue("jekyll/traits/{trait}.md")
system(str_glue("mkdir -p jekyll/traits/{trait}"))

cat("---\n", 'title: "', tbl_traits$NAME[i], '"\npermalink: traits/', trait, "/\nlayout: trait\n", "---\n\n", sep = "", file = fout)
cat("## [Hub]({{ site.baseurl }}) : [Traits]({{ site.baseurl }}traits/)\n\n", file = fout, append = TRUE)
cat("# ", tbl_traits$NAME[i], "\n", sep = "", file = fout, append = TRUE)
cat("`" , length(top_models), " significantly associated models · ", length(top_genes), " unique genes`.\n\n", sep = "", file = fout, append = TRUE)

# ---- Get clumped and conditional loci
cur_clumps <- read_tsv(str_glue("{tbl_traits$OUTPUT[i]}.post.report"), col_types = "ciiiiidddd") |>
    arrange(CHR, P0) |>
    mutate(VAR.EXP = round(VAR.EXP * 100, 0),
           link = str_glue("*[{seq_len(n())}]({{{{ site.baseurl }}}}traits/{trait}/{seq_len(n())})*"),
           genes = "")

# load clumped genes
clump_mod <- vector()
for (ii in seq_len(nrow(cur_clumps))) {
    cur_genes_tbl <- read_tsv(str_glue("{cur_clumps$FILE[ii]}.genes"),
                              col_types = "ccciiidcdcdddiicdddddddddld") |>
        mutate(num = 1:n(),
               gene_id = gene_id(ID),
               link = str_glue("*[{gene_names[gene_id]}]({{{{ site.baseurl }}}}genes/{gene_id})*")) |>
    left_join(select(tbl_panels, PANEL, TISSUE, MODALITY), by = "PANEL")
    
    cur_genes <- sort(unique(cur_genes_tbl$gene_id[cur_genes_tbl$JOINT]))
    cur_genes <- str_c("*", str_c(str_glue("[{gene_names[cur_genes]}]({{{{ site.baseurl }}}}genes/{cur_genes})"), collapse = " "), "*")
    cur_clumps$genes[ii] <- cur_genes
    
    clump_mod <- unique(c(clump_mod, cur_genes_tbl$FILE[cur_genes_tbl$JOINT]))
    
    fout_clump <- str_glue("jekyll/traits/{trait}/{ii}.md")
    cat(str_glue('---\ntitle: "{tbl_traits$NAME[i]}"\npermalink: traits/{trait}/{ii}/\nlayout: locus\n---\n\n'), sep = "", file = fout_clump)
    cat("## [Hub]({{ site.baseurl }}) : [Traits]({{ site.baseurl }}traits) : ", tbl_traits$link[i], " : ", sep = "", file = fout_clump, append = TRUE)
    if (ii > 1) {
        cat(str_glue(" [ ← ]({{{{ site.baseurl }}}}traits/{trait}/{ii-1}) "), sep = "", file = fout_clump, append = TRUE)
    }
    if (ii < nrow(cur_clumps)) {
        cat(str_glue(" [ → ]({{{{ site.baseurl }}}}traits/{trait}/{ii+1})"), sep = "", file = fout_clump, append = TRUE)
    }
    pos0 <- formatC(cur_clumps$P0[ii], format = "f", big.mark = ",", drop0trailing = TRUE)
    pos1 <- formatC(cur_clumps$P1[ii], format = "f", big.mark = ",", drop0trailing = TRUE)
    cat(str_glue("\n\n# chr{cur_clumps$CHR[ii]}:{pos0}-{pos1}\n\n"), sep = "", file = fout_clump, append = TRUE)
    cat("`Best TWAS P=", cur_clumps$BEST.TWAS.P[ii],
        " · Best GWAS P=", cur_clumps$BEST.SNP.P[ii],
        " conditioned to ", cur_clumps$COND.SNP.P[ii],
        "`\n\n", sep = "", file = fout_clump, append = TRUE)
    
    system(str_glue("cp {cur_clumps$FILE[ii]}.cond.csv jekyll/traits/{trait}/{ii}.cond.csv"))
    cat("<script>", "\n", 'Plotly.d3.csv("../', ii, '.cond.csv"', ", function(data){ processData(data) } );", "\n", '</script><div id="graph"></div>', "\n", sep = "", file = fout_clump, append = TRUE)
    # BCAC.1.post.loc_10.cond.csv
    
    cat("\n### Associated models\n\n", sep = "", file = fout_clump, append = TRUE)
    cat("| # | Tissue | Gene | Modality | RNA phenotype | h2 | eQTL R2 | model | # weights | model R2 | model R2 P | eQTL GWAS Z | TWAS Z | TWAS P | Top SNP corr | PP3 | PP4 | joint |", "| --- |", sep = "\n", file = fout_clump, append = TRUE)
    
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
    cat("{: #models}\n\n", file = fout_clump, append = TRUE)
}
cat("\n### Significant Loci\n\n", sep = "", file = fout, append = TRUE)
cat("| # | chr | p0 | p1 | # assoc genes | # joint models | best TWAS P | best SNP P | cond SNP P | % var exp | joint genes |\n| --- |\n", sep = "", file = fout, append = TRUE)
cur_clumps |>
    select(link, CHR, P0, P1, HIT.GENES, JOINT.GENES, BEST.TWAS.P, BEST.SNP.P, COND.SNP.P, VAR.EXP, genes) |>
    as.data.frame() |>
    format(digits = 2) |>
    write.table(quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " | ", file = fout, append = TRUE)
cat("{: #loci}\n\n", file = fout, append = TRUE)

# ---- Get pleiotropic loci
n_pleiot <- 0
for (ii in 1:nrow(tbl_traits)) {
    # cat(ii, "\n")
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
            genes_link <- str_c("*", str_c(str_glue("[{gene_names[genes]}]({{{{ site.baseurl }}}}genes/{genes})"), collapse=" "), "*")
            n_genes_twas <- pleio |>
                filter(TWAS.P.y < 0.05 / n) |>
                with(length(unique(gene_id(ID))))
            df_tmp <- data.frame("link" = tbl_traits$link[ii],
                                 "chisq.ratio" = round(mean(pleio$TWAS.Z.y^2, na.rm = TRUE) / tbl_traits$AVG.CHISQ[ii], 2),
                                 "num.genes" = length(genes),
                                 "num.genes.twas" = n_genes_twas,
                                 "pct.genes.twas" = round(100 * n_genes_twas / tbl_traits$NUM.JOINT.GENES[i], 1), # Replaced 3 with i
                                 "corr" = round(tst$est, 2),
                                 "p.val" = tst$p.value,
                                 "genes" = genes_link)
            if (n_pleiot == 0) {
                df_pleiot <- df_tmp
            } else {
                df_pleiot <- rbind(df_pleiot, df_tmp)
            }
            n_pleiot <- nrow(df_pleiot)
        }
    }
}

cat("### Pleiotropic Associations\n\n", sep = "", file = fout, append = TRUE)
if (n_pleiot != 0) {
    cat("| Trait | chisq ratio | # genes<sup>+</sup> | # genes<sup>++</sup> | % genes<sup>++</sup> | corr | corr P | genes |", "| --- |", sep = "\n", file = fout, append = TRUE)
    df_pleiot |>
        replace_na(list(pct.genes.twas = 0)) |>
        as.data.frame() |>
        format(digits = 2) |>
        write.table(quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " | ", file = fout, append = TRUE)
    cat("{: #pleiotropic}\n\n", file = fout, append = TRUE)
}

cat("### Associations by panel\n\n", sep = "", file = fout, append = TRUE)
cat("| tissue | modality | # hits | % hits/tests | avg chisq |", "| --- |", sep = "\n", file = fout, append = TRUE)
df_cur_models |>
    select(TISSUE, MODALITY, num.hits, pct.hits, avg.chisq) |>
    as.data.frame() |>
    format(digits = 2) |>
    write.table(quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " | ", file = fout, append = TRUE)
cat("{: #panels}\n\n", file = fout, append = TRUE)
