# Load the results of a single trait and generate a markdown file for the trait page and each of its hit loci

# Inputs:
# data/genes.par
# data/traits.par
# data/panels.par
# data/traits.par.nfo
# data/twas_out/{trait}.dat
# data/twas_out/{trait}.dat.post.report
# data/twas_out/{trait}.top

# Outputs:
# jekyll/traits/{trait}.md
# jekyll/traits/{trait}/{locus}.md

suppressPackageStartupMessages(library(tidyverse))

arg <- commandArgs(trailingOnly = TRUE)
i_trait <- as.numeric(arg[1])

# Load gene names since Ensembl IDs were used for TWAS
gene_names <- read_tsv("data/genes.par", col_types = "cc") |>
    deframe()

cat("reading data/traits.par\n")
tbl_traits <- read_tsv("data/traits.par", col_types = "ccicicc")
n_traits <- nrow(tbl_traits)
df_traits <- data.frame("name" = tbl_traits$NAME,
                        "n" = tbl_traits$N,
                        "num.genes" = rep(NA, n_traits),
                        "num.models" = rep(NA, n_traits),
                        "ref" = tbl_traits$REF,
                        "year" = tbl_traits$YEAR,
                        row.names = tbl_traits$ID)
df_traits$link <- str_glue("[{df_traits$name}]({{{{ site.baseurl }}}}traits/{rownames(df_traits)})")

cat("reading data/panels.par\n")
tbl_panels <- read_tsv("data/panels.par", col_types = "ccccci")
n_models <- nrow(tbl_panels)

cat("reading data/traits.par.nfo\n")
traits_nfo <- read_tsv("data/traits.par.nfo", col_types = "ciiid")
m <- match(rownames(df_traits), traits_nfo$ID)
traits_nfo <- traits_nfo[m, ]

# Count the number of significant genes for the trait
i <- i_trait

cat("reading", tbl_traits$OUTPUT[i], "\n")
cur <- read_tsv(
    tbl_traits$OUTPUT[i],
    col_types = cols(PANEL = "c", FILE = "c", ID = "c", TWAS.Z = "d", TWAS.P = "d", .default = "-")
) |>
    filter(!is.na(TWAS.P))

n <- nrow(cur)
top_models <- which(as.numeric(cur$TWAS.P) < 0.05 / n)
n_top <- length(top_models)

top_genes <- unique(cur$ID[top_models])

df_traits$num.models[i] <- length(top_models)
df_traits$num.genes[i] <- length(top_genes)

df_cur_models <- data.frame("tissue" = tbl_panels$TISSUE,
                            "modality" = tbl_panels$MODALITY,
                            "num.hits" = rep(0, n_models),
                            "pct.hits" = rep(0, n_models),
                            "avg.chisq" = rep(NA, n_models),
                            row.names = tbl_panels$ID)
for (m in 1:n_models) {
    keep <- cur$PANEL == tbl_panels$PANEL[m]
    df_cur_models$avg.chisq[m] <- mean(as.numeric(cur$TWAS.Z[keep])^2, na.rm = TRUE)
    df_cur_models$num.hits[m] <- sum(as.numeric(cur$TWAS.P[keep]) < 0.05 / n, na.rm = TRUE)
    df_cur_models$pct.hits[m] <- 100 * mean(as.numeric(cur$TWAS.P[keep]) < 0.05 / n, na.rm = TRUE)
}

# ---- PRINT TRAIT PAGE

fout <- str_glue("jekyll/traits/{tbl_traits$ID[i]}.md")
system(str_glue("mkdir -p jekyll/traits/{tbl_traits$ID[i]}"))

cat("---\n", 'title: "', tbl_traits$NAME[i], '"\npermalink: traits/', tbl_traits$ID[i], "/\nlayout: trait\n", "---\n\n", sep = "", file = fout)
cat("## [Hub]({{ site.baseurl }}) : [Traits]({{ site.baseurl }}traits/)\n\n", file = fout, append = TRUE)
cat("# ", tbl_traits$NAME[i], "\n", sep = "", file = fout, append = TRUE)
cat("`" , df_traits$num.models[i], " significantly associated models · ", df_traits$num.genes[i], " unique genes`.\n\n", sep = "", file = fout, append = TRUE)

# ---- Get clumped and conditional loci
cur_clumps <- read_tsv(str_glue("{tbl_traits$OUTPUT[i]}.post.report"), col_types = "ciiiiidddd")

cat("\n### Significant Loci\n\n", sep = "", file = fout, append = TRUE)
cat("| # | chr | p0 | p1 | # assoc genes | # joint genes | best TWAS P | best SNP P | cond SNP P | % var exp | joint genes |\n| --- |\n", sep = "", file = fout, append = TRUE)
cur_clumps$VAR.EXP <- round(cur_clumps$VAR.EXP * 100, 0)
cur_clumps <- cur_clumps[order(cur_clumps$CHR, cur_clumps$P0), ]
cur_clumps_files <- cur_clumps$FILE
cur_clumps$FILE <- 1:nrow(cur_clumps)
cur_clumps$FILE <- str_glue("*[{cur_clumps$FILE}]({{{{ site.baseurl }}}}traits/{rownames(df_traits)[i]}/{cur_clumps$FILE})*")
cur_clumps$GENES <- rep("", nrow(cur_clumps))

# load clumped genes
clump_mod <- vector()
for (ii in 1:nrow(cur_clumps)) {
    cur_genes_tbl <- read.table(paste(cur_clumps_files[ii], ".genes", sep = ""), head = TRUE, as.is = TRUE)
    cur_genes_tbl$NUM <- 1:nrow(cur_genes_tbl)
    cur_genes_tbl$GENE <- paste("*[", gene_names[cur_genes_tbl$ID], "]({{ site.baseurl }}genes/", cur_genes_tbl$ID, ")*", sep = "")
    
    m <- match(cur_genes_tbl$PANEL, tbl_panels$PANEL)
    cur_genes_tbl$TISSUE <- tbl_panels$TISSUE[m]
    cur_genes_tbl$MODALITY <- tbl_panels$MODALITY[m]
    cur_genes <- sort(unique(cur_genes_tbl$ID[cur_genes_tbl$JOINT]))
    cur_genes <- paste("*", paste("[", gene_names[cur_genes], "]({{ site.baseurl }}genes/", cur_genes, ")", sep = "", collapse = " "), "*", sep = "")
    cur_clumps$GENES[ii] <- cur_genes
    
    clump_mod <- unique(c(clump_mod,cur_genes_tbl$FILE[cur_genes_tbl$JOINT]))
    
    fout_clump <- str_glue("jekyll/traits/{tbl_traits$ID[i]}/{ii}.md")
    cat(str_glue('---\ntitle: "{tbl_traits$NAME[i]}"\npermalink: traits/{tbl_traits$ID[i]}/{ii}/\nlayout: locus\n---\n\n'), sep="", file = fout_clump)
    cat("## [Hub]({{ site.baseurl }}) : [Traits]({{ site.baseurl }}traits) : ", df_traits$link[i], " : ", sep = "", file = fout_clump, append = TRUE)
    if (ii > 1) {
        cat(" [ ← ]({{ site.baseurl }}traits/", tbl_traits$ID[i], "/", (ii-1), ") ", sep = "", file = fout_clump, append = TRUE)
    }
    if (ii < nrow(cur_clumps)) {
        cat(" [ → ]({{ site.baseurl }}traits/", tbl_traits$ID[i], "/", (ii+1), ")", sep = "", file = fout_clump, append = TRUE)
    }
    pos0 <- formatC(cur_clumps$P0[ii], format = "f", big.mark = ",", drop0trailing = TRUE)
    pos1 <- formatC(cur_clumps$P1[ii], format = "f", big.mark = ",", drop0trailing = TRUE)
    cat(str_glue("\n\n# chr{cur_clumps$CHR[ii]}:{pos0}-{pos1}\n\n"), sep = "", file = fout_clump, append = TRUE)
    cat("`Best TWAS P=", cur_clumps$BEST.TWAS.P[ii],
        " · Best GWAS P=", cur_clumps$BEST.SNP.P[ii],
        " conditioned to ", cur_clumps$COND.SNP.P[ii],
        "`\n\n", sep = "", file = fout_clump, append = TRUE)
    
    system(paste("cp ", cur_clumps_files[ii], ".cond.csv jekyll/traits/", tbl_traits$ID[i], "/", ii, ".cond.csv", sep = ""))
    cat("<script>", "\n", 'Plotly.d3.csv("../', ii, '.cond.csv"', ", function(data){ processData(data) } );", "\n", '</script><div id="graph"></div>', "\n", sep = "", file = fout_clump, append = TRUE)
    # BCAC.1.post.loc_10.cond.csv
    
    cat("\n### Associated models\n\n", sep = "", file = fout_clump, append = TRUE)
    cat("| # | Tissue | Modality | Gene | h2 | eQTL R2 | model | # weights | model R2 | model R2 P | eQTL GWAS Z | TWAS Z | TWAS P | Top SNP corr | PP3 | PP4 | joint |", "| --- |", sep = "\n", file = fout_clump, append = TRUE)
    
    cur_genes_tbl$COLOC.PP3 <- round(cur_genes_tbl$COLOC.PP3, 2)
    cur_genes_tbl$COLOC.PP4 <- round(cur_genes_tbl$COLOC.PP4, 2)
    cur_genes_tbl$MODELCV.R2 <- round(cur_genes_tbl$MODELCV.R2, 2)
    cur_genes_tbl$HSQ <- round(cur_genes_tbl$HSQ, 2)
    cur_genes_tbl$EQTL.R2 <- round(cur_genes_tbl$EQTL.R2, 2)
    
    write.table(
        format(cur_genes_tbl[, c("NUM", "TISSUE", "MODALITY", "GENE", "HSQ", "EQTL.R2", "MODEL", "NWGT", "MODELCV.R2", "MODELCV.PV", "EQTL.GWAS.Z", "TWAS.Z", "TWAS.P", "TOP.SNP.COR", "COLOC.PP3", "COLOC.PP4", "JOINT")], digits = 2),
        quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " | ", file = fout_clump, append = TRUE
    )
    cat("{: #models}\n\n", file = fout_clump, append = TRUE)
}
write.table(format(as.data.frame(cur_clumps), digits = 2), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " | ", file = fout, append = TRUE)
cat("{: #loci}\n\n", file = fout, append = TRUE)

# # ---- Get pleiotropic loci
# n.pleiot <- 0
# file.top <- paste("data/twas_out/", tbl_traits$ID[i], ".top", sep = "")
# write.table(cur$FILE[!is.na(match(cur$FILE, clump_mod))], quote = FALSE, row.names = FALSE, col.names = "FILE", file = file.top)
# 
# for (ii in 1:n_traits) {
#     cat(ii, "\n")
#     if (ii != i) {
#         system(paste("cat ", file.top, " | ~/tools/search ", tbl_traits$OUTPUT[ii], " 2 >", file.top, ".tmp", sep = ""))
#         other.cur <- read.table(paste(file.top, ".tmp", sep = ""), as.is = TRUE, head = TRUE, sep = "\t")
#         m <- match(other.cur$FILE, cur$FILE[top_models])
#         other.cur <- other.cur[!is.na(m), ]
#         m <- m[!is.na(m)]
#         other.this <- (cur[top_models, ])[m, ]
#         keep <- which(as.numeric(other.cur$TWAS.P) < 0.05 / n_top)
#         if (length(keep) > 0) {
#             if(length(keep) >= 4) {
#                 tst <- cor.test(as.numeric(other.this$TWAS.Z[keep]), as.numeric(other.cur$TWAS.Z[keep]))
#             } else {
#                 tst <- data.frame("est" = 0, "p.value" = 1)
#             }
#             genes <- sort(unique(other.this$ID[keep]))
#             genes.link <- paste("*", paste("[", gene_names[genes], "]({{ site.baseurl }}genes/", genes, ")", sep = "", collapse=" "), "*", sep = "")
#             num.genes.twas <- length(unique(other.this$ID[which(as.numeric(other.cur$TWAS.P) < 0.05 / n)]))
#             df.tmp <- data.frame("link" = df_traits$link[ii],
#                                  "chisq.ratio" = round(mean(as.numeric(other.cur$TWAS.Z)^2, na.rm = TRUE) / traits_nfo$AVG.CHISQ[ii], 2),
#                                  "num.genes" = length(unique(other.this$ID[keep])),
#                                  "num.genes.twas" = num.genes.twas,
#                                  "pct.genes.twas" = round(100 * num.genes.twas / traits_nfo$NUM.JOINT.GENES[3], 1),
#                                  "corr" = round(tst$est, 2),
#                                  "p.val" = tst$p.value,
#                                  "genes" = genes.link)
#             if (n.pleiot == 0) {
#                 df.pleiot <- df.tmp
#             } else {
#                 df.pleiot <- rbind(df.pleiot, df.tmp)
#             }
#             n.pleiot <- nrow(df.pleiot)
#         }
#     }
# }

# cat("### Pleiotropic Associations\n\n", sep = "", file = fout, append = TRUE)
# if (n.pleiot != 0) {
#     cat("| Trait | chisq ratio | # genes<sup>+</sup> | # genes<sup>++</sup> | % genes<sup>++</sup> | corr | corr P | genes |", "| --- |", sep = "\n", file = fout, append = TRUE)
#     df.pleiot$pct.genes.twas[is.na(df.pleiot$pct.genes.twas)] <- 0
#     write.table(format(df.pleiot, digits = 2), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " | ", file = fout, append = TRUE)
#     cat("{: #pleiotropic}\n\n", file = fout, append = TRUE)
# }
# 
cat("### Associations by panel\n\n", sep = "", file = fout, append = TRUE)
cat("| tissue | modality | # hits | % hits/tests | avg chisq |", "| --- |", sep = "\n", file = fout, append = TRUE)
write.table(format(df_cur_models, digits = 2), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " | ", file = fout, append = TRUE)
cat("{: #panels}\n\n", file = fout, append = TRUE)
