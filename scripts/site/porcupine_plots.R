# Create porcupine plot for each trait's TWAS hits

# Inputs:
# data/traits.par
# data/sequence_report.tsv
# data/all_models.par
# data/twas_out/{trait}.all.tsv
# data/twas_out/{trait}.post.tsv

# Outputs:
# jekyll/assets/images/porcupine/{trait}.png

suppressPackageStartupMessages(library(tidyverse))

process_twas_output <- function(df, chr_len) {
  df |>
    filter(!is.na(TWAS.P)) |>
    # separate_wider_delim(PANEL, ".", names = c("ratgtex", "tissue", "modality")) |>
    mutate(
      # gene_id = str_replace(ID, "__.*", ""),
      chr_num = str_replace(CHR, "chr", "") |> as.integer(),
      logp = -log10(TWAS.P),
      gpos = P0 + cumsum(c(0, chr_len))[chr_num],
      gcolor = as.factor((chr_num - 1) %% 2)
    )
}

process_twas_post <- function(df, chr_len) {
  # Plot midpoint of locus since the full range is typically smaller than a point on these plots
  df |>
    mutate(
      chr_num = str_replace(CHR, "chr", "") |> as.integer(),
      logp = -log10(BEST.TWAS.P),
      gpos = ((P0 + P1) / 2) + cumsum(c(0, chr_len))[chr_num],
      gcolor = as.factor((chr_num - 1) %% 2)
    )
}

porcupine_plot <- function(df_all, df_loci, chr_len, pval_threshold) {
  label_locs <- cumsum(c(0, chr_len[1:19])) + chr_len / 2
  grid_locs <- cumsum(c(0, chr_len))
  points <- process_twas_output(df_all, chr_len)
  loci <- process_twas_post(df_loci, chr_len)
  ggplot(points, aes(x = gpos, y = logp, color = gcolor)) +
    geom_hline(yintercept = -log10(pval_threshold), color = "red", lty = "22") +
    geom_point(size = 0.5, show.legend = FALSE) +
    geom_point(aes(fill = HIT.GENES), color = "black", shape = 21, data = loci, size = 2) +
    expand_limits(x = c(0, sum(chr_len))) +
    scale_x_continuous(breaks = label_locs, labels = 1:20, expand = c(0.02, 0),
                       minor_breaks = grid_locs)  +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    scale_color_manual(values = c("#555555", "#aaaaaa")) +
    scale_fill_viridis_c(trans = "log", breaks = c(1, 3, 10, 30, 100, 300, 1000)) +
    theme_classic() +
    theme(
      legend.key.width = unit(10, "pt"),
    ) +
    xlab("Chromosome") +
    ylab(expression(-log[10]*P)) +
    labs(fill = "Signif.\ngenes\nin locus")
}

traits <- read_tsv("data/traits.par", col_types = "-c-----") |> pull()

chrom_lengths <- read_tsv(
  "data/sequence_report.tsv",
  col_types = cols(`Seq length` = "i", `Sequence name` = "c")
) |>
  select(chrom = `Sequence name`, length = `Seq length`) |>
  filter(chrom %in% str_c("chr", 1:20)) |>
  mutate(chrom = str_replace(chrom, "chr", "") |> as.integer()) |>
  deframe()

# This is similar to row count of TWAS output, but some phenotypes are missing
# from output, and p-value threshold is calculated from all models.
n_models <- length(read_lines("data/all_models.par")) - 1

for (trait in traits) {
  df_all <- read_tsv(
    str_glue("data/twas_out/{trait}.all.tsv"),
    col_types = cols(PANEL = "c", ID = "c", CHR = "c", P0 = "i", TWAS.P = "d", .default = "-")
  )
  df_loci <- read_tsv(
    str_glue("data/twas_out/{trait}.post.tsv"),
    col_types = cols(CHR = "c", P0 = "i", P1 = "i", HIT.GENES = "i", BEST.TWAS.P = "d", .default = "-")
  )
  p <- porcupine_plot(df_all, df_loci, chrom_lengths, 0.05 / n_models)
  ggsave(str_glue("jekyll/assets/images/porcupine/{trait}.png"), p,
         width = 8, height = 3, device = png)
}
