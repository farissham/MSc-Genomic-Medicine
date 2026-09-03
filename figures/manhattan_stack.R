suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
})

base <- "/Users/USER-PC/Desktop/hcm_dcm_output"
hub_dir <- file.path(base, "results_comparison_all4/harmonise")

datasets <- list(
  dcm        = list(panel_letter = "d",
                    colors = c("#e87ba4", "#F6CADB"),
                    hub = file.path(hub_dir, "zheng_dcm.harmonised.tsv.gz"),
                    magma = file.path(base, "results_zheng/magma/zheng_dcm.genes.out"),
                    susie_dir = file.path(base, "results_zheng/susie"),         prefix = "zheng_dcm"),
  hcm        = list(panel_letter = "c",
                    colors = c("#4a3aa7", "#B7B0DC"),
                    hub = file.path(hub_dir, "tadros_hcm.harmonised.tsv.gz"),
                    magma = file.path(base, "results_tadros/magma/tadros_hcm.genes.out"),
                    susie_dir = file.path(base, "results_tadros/susie"),        prefix = "tadros_hcm"),
  kramarenko = list(panel_letter = "b",
                    colors = c("#005F73", "#99BFC7"),
                    hub = file.path(hub_dir, "kramarenko_ccmtag.harmonised.tsv.gz"),
                    magma = file.path(base, "results_kramarenko/magma/kramarenko_ccmtag.genes.out"),
                    susie_dir = file.path(base, "results_kramarenko/susie"),    prefix = "kramarenko_ccmtag"),
  primary    = list(panel_letter = "a",
                    colors = c("#2a78d6", "#a9c8ef"),
                    hub = file.path(hub_dir, "hcm_dcm.harmonised.tsv.gz"),
                    magma = file.path(base, "results_hcm_dcm/magma/hcm_dcm.genes.out"),
                    susie_dir = file.path(base, "results_hcm_dcm/susie"),       prefix = "hcm_dcm")
)

chr_order  <- as.character(1:22)
genome_wide_sig <- 5e-8

for (key in names(datasets)) {

  ds <- datasets[[key]]
  cat("\n==", ds$label, "==\n")

  hub <- fread(ds$hub, select = c("chr", "pos", "p"))
  hub[, chr := as.character(chr)]
  hub <- hub[chr %in% chr_order]
  hub[, chr := factor(chr, levels = chr_order)]

  chr_len <- hub[, .(chr_max = max(pos)), by = chr][order(chr)]
  chr_len[, chr_add := cumsum(as.numeric(shift(chr_max, fill = 0)))]
  hub <- merge(hub, chr_len, by = "chr")
  hub[, pos_cum := pos + chr_add]
  hub[, negLog10P := -log10(p)]
  hub[, col_group := as.integer(chr) %% 2]

  axis_df <- hub[, .(center = mean(pos_cum)), by = chr][order(chr)]
  n_gw_hits <- sum(hub$p < genome_wide_sig)
  cat("  SNPs plotted:", nrow(hub), "| genome-wide p <", genome_wide_sig, "hits:", n_gw_hits, "\n")

  p_manhattan <- ggplot(hub, aes(x = pos_cum, y = negLog10P, color = factor(col_group))) +
    geom_point(size = 0.6, alpha = 0.6) +
    geom_hline(yintercept = -log10(genome_wide_sig), linetype = "dashed", color = "#e34948", linewidth = 0.4) +
    scale_color_manual(values = ds$colors, guide = "none")  +
    scale_x_continuous(labels = axis_df$chr, breaks = axis_df$center, expand = c(0.01, 0.01),
                        limits = c(0, max(hub$pos_cum))) +
    coord_cartesian(ylim = c(0, 120)) +
    labs(x = NULL, y = expression(-log[10](italic(P))), title = "GWAS Manhattan (per-SNP)") +
    theme_classic(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(size = 10, face = "bold"))

  cs_files <- list.files(ds$susie_dir, pattern = "\\.credible_sets\\.csv$", full.names = TRUE)
  susie_list <- list()

  for (f in cs_files) {
    lines <- readLines(f, n = 2)
    if (length(lines) < 2 || grepl("no credible sets found", lines[2])) next

    fname <- basename(f)
    m <- regmatches(fname, regexec(paste0("^", ds$prefix, "_([0-9]+)_([0-9]+)_([0-9]+)\\.credible_sets\\.csv$"), fname))[[1]]
    if (length(m) == 0) next
    locus_chr <- m[2]

    d <- fread(f)
    for (i in seq_len(nrow(d))) {
      pos_vals <- as.numeric(strsplit(as.character(d$pos[i]), ",")[[1]])
      pip_vals <- as.numeric(strsplit(as.character(d$pip[i]), ",")[[1]])
      susie_list[[length(susie_list) + 1]] <- data.table(chr = locus_chr, pos = pos_vals, PIP = pip_vals, cs = d$cs[i])
    }
  }

  s <- rbindlist(susie_list)
  s <- s[chr %in% chr_order]
  s[, chr := factor(chr, levels = chr_order)]
  s <- merge(s, chr_len[, .(chr, chr_add)], by = "chr")
  s[, pos_cum := pos + chr_add]
  cat("  SuSiE credible-set variants plotted:", nrow(s), "across", length(unique(paste(s$chr, s$cs))), "credible sets\n")

  p_susie <- ggplot(s, aes(x = pos_cum, y = PIP, color = PIP)) +
    geom_point(size = 1.3, alpha = 0.85) +
    scale_color_gradient(low = "#a9c8ef", high = "#08306b", limits = c(0, 1), name = "PIP") +
    scale_x_continuous(labels = axis_df$chr, breaks = axis_df$center, expand = c(0.01, 0.01),
                       limits = c(0, max(hub$pos_cum))) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = NULL, y = "SuSiE PIP", title = "SuSiE fine-mapping (credible-set variants)") +
    theme_classic(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(size = 10, face = "bold"))

  g <- fread(ds$magma)
  g[, CHR := as.character(CHR)]
  g <- g[CHR %in% chr_order]
  g[, CHR := factor(CHR, levels = chr_order)]
  g[, POS := (START + STOP) / 2]
  g[, negLog10P := -log10(P)]
  g <- merge(g, chr_len[, .(chr, chr_add)], by.x = "CHR", by.y = "chr")
  g[, POS_cum := POS + chr_add]
  g[, col_group := as.integer(CHR) %% 2]

  bonf   <- 0.05 / nrow(g)
  n_hits <- sum(g$P < bonf)
  cat("  MAGMA genes tested:", nrow(g), "| Bonferroni p <", signif(bonf, 3), "| hits:", n_hits, "\n")

  p_magma <- ggplot(g, aes(x = POS_cum, y = negLog10P, color = factor(col_group))) +
    geom_point(size = 0.9, alpha = 0.75) +
    geom_hline(yintercept = -log10(bonf), linetype = "dashed", color = "#e34948", linewidth = 0.4) +
    scale_color_manual(values = ds$colors, guide = "none") +
    scale_x_continuous(labels = axis_df$chr, breaks = axis_df$center, expand = c(0.01, 0.01),
                        limits = c(0, max(hub$pos_cum))) +
    coord_cartesian(ylim = c(0, 15)) +
    labs(x = "Chromosome", y = expression(-log[10](italic(P))), title = "MAGMA gene-based test") +
    theme_classic(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text.x = element_text(size = 7),
          plot.title = element_text(size = 10, face = "bold"))

  if (n_hits > 0) {
    top_hits <- g[P < bonf][order(P)][1:min(.N, 15)]
    p_magma <- p_magma + geom_text_repel(data = top_hits, aes(label = GENE), color = "black",
                                          size = 2.4, max.overlaps = 20, segment.size = 0.2)
  }

  combined <- (p_manhattan / p_susie / p_magma) +
    plot_layout(heights = c(1, 1, 1), guides = "collect") &
    theme(legend.position = "right")
  
  combined <- combined +
    plot_annotation(title = paste0(ds$panel_letter, "  ", ds$label), theme = theme(plot.title = element_text(size = 13, face = "bold")))

  out_png <- file.path(base, paste0("manhattan_stack_", key, ".png"))
  ggsave(out_png, combined, width = 10, height = 10, dpi = 300, bg = "white")
  cat("  Saved:", out_png, "\n")
}
