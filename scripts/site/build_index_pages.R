# Load the data and generate the index pages for the traits, models, genes, and projects.

# Inputs:
# data/traits.par
# data/traits.par.nfo
# data/all_models.par
# data/genes_n_assoc.nfo
# data/genes_n_models.nfo
# data/gene_names.tsv

# Outputs:
# jekyll/traits.md
# jekyll/genes.json
# jekyll/_data/genes_stats.yml

suppressPackageStartupMessages(library(tidyverse))
library(yaml)

gene_id <- function(pheno_ID) {
    # Extract gene ID from RNA phenotype ID
    str_extract(pheno_ID, "^[^:.]+")
}

# ---- PRINT TRAIT INDEX
traits_nfo <- read_tsv("data/traits.par.nfo", col_types = "ciiid")

tbl_traits <- read_tsv("data/traits.par", col_types = "ccicccc") |>
    mutate(
        link = str_glue("[{NAME}]({{{{ site.baseurl }}}}traits/{ID})"),
        data_url = str_glue("{{{{ site.baseurl }}}}data/{ID}.tar.bz2"),
        data_link = str_glue('[ <i class="far fa-file-archive" aria-hidden="true"></i> ]({data_url})'),
        project_link = str_glue("[{PROJECT}]({{{{ site.baseurl }}}}projects/)"),
    ) |>
    left_join(traits_nfo, by = "ID")

fout <- "jekyll/traits.md"
cat("---", "title: Traits", "permalink: traits/", "layout: traits", "---\n", sep = "\n", file = fout)
n_loci <- formatC(sum(tbl_traits$NUM.LOCI), format = "f", big.mark = ",", drop0trailing = TRUE)
n_genes <- formatC(sum(tbl_traits$NUM.GENES), format = "f", big.mark = ",", drop0trailing = TRUE)
cat("{: .text-center }\n### **",
    nrow(tbl_traits), "** traits &middot; **",
    n_loci, "** associated loci &middot; **",
    n_genes, "**  gene/trait associations\n\n",
    sep = "", file = fout, append = TRUE)
cat("| Trait | N | # loci | # indep genes | # total genes | Project | data | ", "| --- |", sep = "\n", file = fout, append = TRUE)
tbl_traits |>
    select(link, N, NUM.LOCI, NUM.JOINT.GENES, NUM.GENES, project_link, data_link) |>
    write.table(quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " | ", file = fout, append = TRUE)

## Make genes.json
all_genes <- read_tsv("data/all_models.par", col_types = cols(ID = "c", .default = "-")) |>
    mutate(gene_id = gene_id(ID)) |>
    with(sort(unique(gene_id)))
genes_n_assoc <- read_delim("data/genes_n_assoc.nfo", delim = " ", col_types = "ci", col_names = FALSE) |>
    deframe()
genes_n_models <- read_delim("data/genes_n_models.nfo", delim = " ", col_types = "ci", col_names = FALSE) |>
    deframe()
# Load gene names since Ensembl IDs were used for TWAS
gene_names <- read_tsv("data/gene_names.tsv", col_types = "cc") |>
    deframe()

df_genes <- tibble(gene = all_genes) |>
    mutate(n.assoc = genes_n_assoc[gene],
           n.models = genes_n_models[gene],
           gene_link = str_glue('<a href=\\"./{gene}\\">{gene_names[gene]}</a>')) |>
    replace_na(list(n.assoc = 0, n.models = 0))

cat('{\n"data":[\n',
    str_c(str_glue('["{df_genes$gene_link}","{df_genes$gene}",{df_genes$n.assoc},{df_genes$n.models}]'), collapse = ",\n"),
    "]\n}",
    sep = "",
    file = "jekyll/genes.json")

# ---- STATS FOR GENE INDEX
gene_stats <- list(
  n_genes = formatC(nrow(df_genes), format = "f", big.mark = ",", drop0trailing = TRUE),
  n_models = formatC(sum(df_genes$n.models), format = "f", big.mark = ",", drop0trailing = TRUE)
)
dir.create("jekyll/_data", showWarnings = FALSE, recursive = TRUE)
yaml::write_yaml(gene_stats, file = "jekyll/_data/genes_stats.yml")
