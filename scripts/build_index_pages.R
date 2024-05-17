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

# Load gene names since Ensembl IDs were used for TWAS
gene_names <- read.table("data/genes.par", as.is = TRUE, head = TRUE, sep = '\t')
gene_names <- setNames(gene_names$NAME, gene_names$ID)

tbl.panels <- read.table("data/panels.par", as.is = TRUE, head = TRUE, sep = '\t')
N.models <- nrow(tbl.panels)
# df.panels = data.frame( "n"=tbl.panels$N , "study"=tbl.panels$STUDY , "tissue"=tbl.panels$TISSUE , "num.hits" = rep(0,N.models) , "num.traits" = rep(0,N.models) , row.names=tbl.panels$ID )
df.panels <- data.frame("n" = tbl.panels$N,
                        "tissue" = tbl.panels$TISSUE,
                        "modality" = tbl.panels$MODALITY,
                        "num.hits" = rep(0, N.models),
                        "num.traits" = rep(0, N.models),
                        row.names = tbl.panels$ID)

tbl.models.pos <- read.table("data/all.models.par", as.is = TRUE, head = TRUE)
m <- match(tbl.models.pos$PANEL, tbl.panels$PANEL)
# tbl.models.pos$PANEL = paste( tbl.panels$STUDY , " | " , tbl.panels$TISSUE , sep='' )[ m ] 
tbl.models.pos$PANEL <- paste(tbl.panels$TISSUE, " | ", tbl.panels$MODALITY, sep = '')[m]

# ---- PRINT TRAIT INDEX
tbl.traits <- read.table("data/traits.par", as.is = TRUE, head = TRUE, sep = '\t', quote = "")
N.traits <- nrow(tbl.traits)
df.traits <- data.frame("type" = tbl.traits$TYPE,
                        "name" = tbl.traits$NAME,
                        "n" = tbl.traits$N,
                        "num.loci" = rep(NA, N.traits),
                        "num.joint.genes" = rep(NA, N.traits),
                        "num.total.genes" = rep(NA, N.traits),
                        "ref" = tbl.traits$REF,
                        "year" = tbl.traits$YEAR,
                        row.names = tbl.traits$ID)
df.traits$link <- paste("[", df.traits$name, "]({{ site.baseurl }}traits/", rownames(df.traits), ")", sep = '')

df.traits$data <- paste("[ <i class=\"far fa-file-archive\" aria-hidden=\"true\"></i> ]({{ site.baseurl }}data/", gsub(".dat", "", gsub("^.+/", "", tbl.traits$OUTPUT)), ".tar.bz2)" , sep = '')

traits.nfo <- read.table("data/traits.par.nfo", as.is = TRUE, head = TRUE)
m <- match(rownames(df.traits), traits.nfo$ID)
traits.nfo <- traits.nfo[m, ]
df.traits$num.loci <- traits.nfo$NUM.LOCI
df.traits$num.joint.genes <- traits.nfo$NUM.JOINT.GENES
df.traits$num.total.genes <- traits.nfo$NUM.GENES

fout <- "jekyll/traits.md"
cat("---", "title: Traits", "permalink: traits/", "layout: traits", "---\n", sep = '\n', file = fout)
cat("# *",
    nrow(df.traits),
    "* traits &middot; *",
    formatC(sum(df.traits$num.loci), format = "f", big.mark = ",", drop0trailing = TRUE),
    "* associated loci &middot; *",
    formatC(sum(df.traits$num.total.genes), format = "f", big.mark = ",", drop0trailing = TRUE),
    "*  gene/trait associations\n\n",
    sep = '', file = fout, append = TRUE)
cat("| Type | Trait | N | # loci | # indep genes | # total genes | Ref. | Year | data | ", "| --- |", sep = '\n', file = fout, append = TRUE)
write.table(
    df.traits[,c("type","link","n","num.loci","num.joint.genes","num.total.genes","ref","year","data")],
    quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ' | ', file = fout, append = TRUE
)

fout <- "jekyll/models.md"
cat("---", "title: Models", "permalink: models/", "layout: about", "---\n", sep = '\n', file = fout)
cat("# Models \n\n", sep = '', file = fout, append = TRUE)

# cat( "| Study | Tissue | N |","| --- |",sep='\n',file=fout,append=T)
# write.table(df.panels[,c("study","tissue","n")],quote=F,row.names=F,col.names=F,sep=' | ',file=fout,append=T)
cat("| Tissue | Modality | N |", "| --- |", sep = '\n', file = fout, append = TRUE)
write.table(df.panels[, c("tissue", "modality", "n")], quote = FALSE,
            row.names = FALSE, col.names = FALSE, sep = ' | ', file = fout, append = TRUE)
# ----

## Make genes.json
uni.genes <- sort(unique(tbl.models.pos$ID))
df.genes <- data.frame("gene" = uni.genes,
                       "n.models" = rep(0, length(uni.genes)),
                       "n.assoc" = rep(0, length(uni.genes)))
# df.genes$link = paste( '*' , paste( "[" , gene_names[df.genes$gene] , "]({{ site.baseurl }}genes/" , df.genes$gene , ")" , sep='' ) , '*' , sep='' )

genes.nfo <- read.table("data/genes.nfo", as.is = TRUE)
df.genes$n.assoc <- genes.nfo[match(df.genes$gene, genes.nfo[, 1]), 2]
genes.nfo <- read.table("data/genes.models.nfo", as.is = TRUE)
df.genes$n.models <- genes.nfo[match(df.genes$gene, genes.nfo[, 1]), 2]
df.genes$n.assoc[is.na(df.genes$n.assoc)] <- 0
df.genes$n.models[is.na(df.genes$n.models)] <- 0

df.json <- df.genes[, c("gene", "n.assoc", "n.models")]
df.json$gene_link <- paste("<em><a href=\\\"./", df.json$gene, "\\\">", gene_names[df.genes$gene], "</a></em>", sep = '')
cat("{\n\"data\":[\n",
    paste("[\"", df.json$gene_link, "\",\"", df.json$gene, "\",", df.json$n.assoc, ",", df.json$n.models, "]", sep = '', collapse = ",\n"),
    "]\n}",
    sep = '',
    file = "jekyll/genes.json")

# ---- PRINT GENE INDEX
fout <- "jekyll/genes.md"
cat("---", "title: Genes", "permalink: genes/", "layout: genes", "---\n", sep = '\n', file = fout)
cat("# *",
    formatC(nrow(df.genes), format = "f", big.mark = ",", drop0trailing = TRUE),
    "* genes &middot; *",
    formatC(sum(df.genes$n.models), format = "f", big.mark = ",", drop0trailing = TRUE),
    "* models\n\n",
    sep = '', file = fout, append = TRUE)
cat("| Gene | ID | # associations | # models |\n", "| --- |\n| |\n", sep = '', file = fout, append = TRUE)
## Table rows get loaded from genes.json instead
#write.table(df.genes[,c("link","n.assoc","n.models")],quote=F,row.names=F,col.names=F,sep=' | ',file=fout,append=T)
cat('{: #genes}\n', file = fout, append = TRUE)
# ----
