# Generate gene pages

# Inputs:
# data/gene_names.tsv
# data/panels.par
# data/all_models.par
# data/traits.par
# data/traits.par.nfo

# Outputs:
# jekyll/genes/{gene}.md

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

# Load gene names since Ensembl IDs were used for TWAS
gene_names <- read_tsv("data/gene_names.tsv", col_types = "cc") |>
    deframe()

tbl_panels <- read_tsv("data/panels.par", col_types = "ccccci")

tbl_models <- read_tsv("data/all_models.par", col_types = "ccciiiiddddddddddddd") |>
    left_join(select(tbl_panels, PANEL, TISSUE, MODALITY), by = "PANEL") |>
    mutate(gene_id = gene_id(ID),
           ID = shorten_ids(ID, MODALITY))

# ---- PRINT TRAIT INDEX

traits_nfo <- read_tsv("data/traits.par.nfo", col_types = "ciiid")

tbl_traits <- read_tsv("data/traits.par", col_types = "ccicicc") |>
    mutate(link = str_glue("[{NAME}]({{{{ site.baseurl }}}}traits/{ID})")) |>
    left_join(traits_nfo, by = "ID")

# iterate and make each gene file
all_genes <- sort(unique(tbl_models$gene_id))

df_genes <- tibble(gene = all_genes) |>
    mutate(n.models = 0,
           link = str_glue("*[{gene}]({{{{ site.baseurl }}}}genes/{gene})*"))

for (i in 1:length(all_genes)) {
    fout <- str_glue("jekyll/genes/{all_genes[i]}.md")
    cat("---\n", "title: ", gene_names[all_genes[i]], "\npermalink: genes/", all_genes[i], "/ \nlayout: gene\n", "---\n\n", sep = "", file = fout)
    cat("## [Hub]({{ site.baseurl }}) : [Genes]({{ site.baseurl }}genes)\n\n", file = fout, append = TRUE)
    cat("# ", gene_names[all_genes[i]], "\n", sep = "", file = fout, append = TRUE)
    
    # WGT	ID	CHR	P0	P1	ID	NSNPS	HSQ	HSQ.SE	HSQ.PV	TOP1.R2	BLUP.R2	ENET.R2	BSLMM.R2	LASSO.R2	TOP1.PV	BLUP.PV	ENET.PV	BSLMM.PV	LASSO.PV
    cat("\n### Models\n\n", sep = "", file = fout, append = TRUE)
    cat("| # | tissue | modality | RNA phenotype | h<sup>2</sup> | h<sup>2</sup> se | h<sup>2</sup> P | eQTL R<sup>2</sup> | BLUP R<sup>2</sup> | ENET R<sup>2</sup> | LASSO R<sup>2</sup> | eQTL P | BLUP P | ENET P | LASSO P |\n| --- |\n", sep = "", file = fout, append = TRUE)
    cur_models <- tbl_models |>
        filter(gene_id == all_genes[i]) |>
        mutate(NUM = 1:n())
    df_genes$n.models[i] <- df_genes$n.models[i] + nrow(cur_models)
    cur_models |>
        select(NUM, TISSUE, MODALITY, ID, HSQ, HSQ.SE, HSQ.PV, TOP1.R2, BLUP.R2,
               ENET.R2, LASSO.R2, TOP1.PV, BLUP.PV, ENET.PV, LASSO.PV) |>
        as.data.frame() |>
        format(digits = 3) |>
        write.table(quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " | ", file = fout, append = TRUE)
    cat("{: #models}\n\n", file = fout, append = TRUE)

    cat("\n### Trait associations\n\n| Trait | Avg chi<sup>2</sup> ratio | Avg chi<sup>2</sup> | Max chi<sup>2</sup> | ",
        str_c(1:df_genes$n.models[i], collapse = " | "),
        " | \n| --- |\n",
        sep = "",
        file = str_glue("jekyll/genes/{all_genes[i]}.md"),
        append = TRUE)
    
    if (i %% 100 == 0) cat(i, "\n")
}

# iterate over each trait and count the number of significant genes
for (i in 1:nrow(tbl_traits)) {
    system(str_glue('bash scripts/add_trait_to_gene_pages.sh {tbl_traits$OUTPUT[i]} "{tbl_traits$link[i]}" {tbl_traits$AVG.CHISQ[i]}\n'))
    cat(i, "\n")
}

for (i in 1:length(all_genes)) {
    cat("{: #assoc}\n", file = str_glue("jekyll/genes/{all_genes[i]}.md"), append = TRUE)
}
