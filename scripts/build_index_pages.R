# Load the data and generate the index pages for the traits, models, and genes.

# Inputs:
# data/genes.par
# data/panels.par
# data/all.models.par
# data/traits.par
# data/traits.par.nfo
# data/genes.nfo
# data/genes.models.nfo

# Outputs:
# jekyll/traits.md
# jekyll/models.md
# jekyll/genes.json
# jekyll/genes.md

suppressPackageStartupMessages(library(tidyverse))

# Load gene names since Ensembl IDs were used for TWAS
gene_names <- read_tsv("data/genes.par", col_types = "cc") |>
    deframe()

tbl_panels <- read_tsv("data/panels.par", col_types = "ccccci")

tbl_models_pos <- read.table("data/all.models.par", as.is = TRUE, head = TRUE)
m <- match(tbl_models_pos$PANEL, tbl_panels$PANEL)
tbl_models_pos$PANEL <- paste(tbl_panels$TISSUE, " | ", tbl_panels$MODALITY, sep = "")[m]

# ---- PRINT TRAIT INDEX
traits_nfo <- read_tsv("data/traits.par.nfo", col_types = "ciiid")

tbl_traits <- read_tsv("data/traits.par", col_types = "ccicicc") |>
    mutate(link = str_glue("[{NAME}]({{{{ site.baseurl }}}}traits/{ID})"),
           data_url = str_glue("{{{{ site.baseurl }}}}data/{ID}.tar.bz2"),
           data_link = str_glue('[ <i class="far fa-file-archive" aria-hidden="true"></i> ]({data_url})')) |>
    left_join(traits_nfo, by = "ID")

fout <- "jekyll/traits.md"
cat("---", "title: Traits", "permalink: traits/", "layout: traits", "---\n", sep = "\n", file = fout)
n_loci <- formatC(sum(tbl_traits$NUM.LOCI), format = "f", big.mark = ",", drop0trailing = TRUE)
n_genes <- formatC(sum(tbl_traits$NUM.GENES), format = "f", big.mark = ",", drop0trailing = TRUE)
cat("# *", nrow(tbl_traits), "* traits &middot; *",
    n_loci, "* associated loci &middot; *",
    n_genes, "*  gene/trait associations\n\n",
    sep = "", file = fout, append = TRUE)
cat("| Type | Trait | N | # loci | # indep genes | # total genes | Ref. | Year | data | ", "| --- |", sep = "\n", file = fout, append = TRUE)
tbl_traits |>
    select(TYPE, link, N, NUM.LOCI, NUM.JOINT.GENES, NUM.GENES, REF, YEAR, data_link) |>
    write.table(quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " | ", file = fout, append = TRUE)

fout <- "jekyll/models.md"
cat("---", "title: Models", "permalink: models/", "layout: about", "---\n", sep = "\n", file = fout)
cat("# Models \n\n", sep = "", file = fout, append = TRUE)

cat("| Tissue | Modality | N |", "| --- |", sep = "\n", file = fout, append = TRUE)
tbl_panels |>
    select(TISSUE, MODALITY, N) |>
    write.table(quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " | ", file = fout, append = TRUE)
# ----

## Make genes.json
uni_genes <- sort(unique(tbl_models_pos$ID))

genes_n_assoc <- read_delim("data/genes.nfo", delim = " ", col_types = "ci", col_names = FALSE) |>
    deframe()
genes_n_models <- read_delim("data/genes.models.nfo", delim = " ", col_types = "ci", col_names = FALSE) |>
    deframe()

df_genes <- tibble(gene = uni_genes) |>#, n.models = 0, n.assoc = 0)
    mutate(n.assoc = genes_n_assoc[gene],
           n.models = genes_n_models[gene],
           gene_link = str_glue('<em><a href=\\"./{gene}\\">{gene_names[gene]}</a></em>')) |>
    replace_na(list(n.assoc = 0, n.models = 0))

cat('{\n"data":[\n',
    str_c(str_glue('["{df_genes$gene_link}","{df_genes$gene}",{df_genes$n.assoc},{df_genes$n.models}]'), collapse = ",\n"),
    "]\n}",
    sep = "",
    file = "jekyll/genes.json")

# ---- PRINT GENE INDEX
fout <- "jekyll/genes.md"
cat("---", "title: Genes", "permalink: genes/", "layout: genes", "---\n", sep = "\n", file = fout)
n_genes <- formatC(nrow(df_genes), format = "f", big.mark = ",", drop0trailing = TRUE)
n_models <- formatC(sum(df_genes$n.models), format = "f", big.mark = ",", drop0trailing = TRUE)
cat(str_glue("# *{n_genes}* genes &middot; *{n_models}* models\n\n"), sep = "", file = fout, append = TRUE)
cat("| Gene | ID | # associations | # models |\n", "| --- |\n| |\n", sep = "", file = fout, append = TRUE)
## Table rows get loaded from genes.json instead
#write.table(df_genes[,c("link","n.assoc","n.models")],quote=F,row.names=F,col.names=F,sep=' | ',file=fout,append=T)
cat("{: #genes}\n", file = fout, append = TRUE)
# ----
