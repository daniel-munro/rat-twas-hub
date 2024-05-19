# Generate gene pages

# Inputs:
# data/genes.par
# data/panels.par
# data/all.models.par
# data/traits.par
# data/traits.par.nfo

# Outputs:
# jekyll/genes/{gene}.md

suppressPackageStartupMessages(library(tidyverse))

# Load gene names since Ensembl IDs were used for TWAS
gene_names <- read_tsv("data/genes.par", col_types = "cc") |>
    deframe()

tbl_panels <- read_tsv("data/panels.par", col_types = "ccccci")

tbl_models <- read_tsv("data/all.models.par", col_types = "ccciiiiddddddddddddd") |>
    left_join(select(tbl_panels, PANEL, TISSUE, MODALITY), by = "PANEL")

# ---- PRINT TRAIT INDEX

traits_nfo <- read_tsv("data/traits.par.nfo", col_types = "ciiid")

tbl_traits <- read_tsv("data/traits.par", col_types = "ccicicc") |>
    mutate(link = str_glue("[{NAME}]({{{{ site.baseurl }}}}traits/{ID})")) |>
    left_join(traits_nfo, by = "ID")

# iterate and make each gene file
uni_genes <- sort(unique(tbl_models$ID))

df_genes <- tibble(gene = uni_genes) |>
    mutate(n.models = 0,
           link = str_glue("*[{gene}]({{{{ site.baseurl }}}}genes/{gene})*"))

for (i in 1:length(uni_genes)) {
    fout <- str_glue("jekyll/genes/{uni_genes[i]}.md")
    cat("---\n", "title: ", gene_names[uni_genes[i]], "\npermalink: genes/", uni_genes[i], "/ \nlayout: gene\n", "---\n\n", sep = "", file = fout)
    cat("## [Hub]({{ site.baseurl }}) : [Genes]({{ site.baseurl }}genes)\n\n", file = fout, append = TRUE)
    cat("# ", gene_names[uni_genes[i]], "\n", sep = "", file = fout, append = TRUE)
    
    # WGT	ID	CHR	P0	P1	ID	NSNPS	HSQ	HSQ.SE	HSQ.PV	TOP1.R2	BLUP.R2	ENET.R2	BSLMM.R2	LASSO.R2	TOP1.PV	BLUP.PV	ENET.PV	BSLMM.PV	LASSO.PV
    cat("\n### Models\n\n", sep = "", file = fout, append = TRUE)
    cat("| # | tissue | modality | h2 | h2 se | h2 P | eQTL R2 | BLUP R2 | ENET R2 | BSLMM R2 | LASSO R2 | eQTL P | BLUP P | ENET P | BSLMM P | LASSO P |\n| --- |\n", sep = "", file = fout, append = TRUE)
    cur_models <- tbl_models |>
        filter(ID == uni_genes[i]) |>
        mutate(NUM = 1:n())
    df_genes$n.models[i] <- df_genes$n.models[i] + nrow(cur_models)
    cur_models |>
        select(NUM, TISSUE, MODALITY, HSQ, HSQ.SE, HSQ.PV, TOP1.R2, BLUP.R2,
               ENET.R2, BSLMM.R2, LASSO.R2, TOP1.PV, BLUP.PV, ENET.PV, BSLMM.PV,
               LASSO.PV) |>
        as.data.frame() |>
        format(digits = 3) |>
        write.table(quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " | ", file = fout, append = TRUE)
    cat("{: #models}\n\n", file = fout, append = TRUE)

    cat("\n### Trait associations\n\n | Trait | Avg chi2 ratio | Avg chi2 | Max chi2 | ",
        str_c(1:df_genes$n.models[i], collapse = " | "),
        "| \n | --- |\n",
        sep = "",
        file = str_glue("jekyll/genes/{uni_genes[i]}.md"),
        append = TRUE)
    
    if (i %% 100 == 0) cat(i, "\n")
}

# iterate over each trait and count the number of significant genes
for (i in 1:nrow(tbl_traits)) {
    system(str_glue('bash scripts/add_trait_to_gene_pages.sh {tbl_traits$OUTPUT[i]} "{tbl_traits$link[i]}" {tbl_traits$AVG.CHISQ[i]}\n'))
    cat(i, "\n")
}

for (i in 1:length(uni_genes)) {
    cat("{: #assoc}\n", file = str_glue("jekyll/genes/{uni_genes[i]}.md"), append = TRUE)
}
