suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

base <- "/Users/USER-PC/Desktop/hcm_dcm_output"

datasets <- list(
  primary    = list(label = "Primary",    dir = "results_hcm_dcm",    color = "#2a78d6"),
  kramarenko = list(label = "Reconstructed", dir = "results_kramarenko", color = "#005F73"),
  hcm        = list(label = "HCM",        dir = "results_tadros",     color = "#4a3aa7"),
  dcm        = list(label = "DCM",        dir = "results_zheng",      color = "#e87ba4")
)

read_coloc_dir <- function(dirname) {
  files <- list.files(file.path(base, dirname, "coloc"), pattern = "\\.coloc_summary\\.tsv$", full.names = TRUE)
  rbindlist(lapply(files, fread), fill = TRUE)
}

make_quadrant_plot <- function(key, ref_suffix, ref_label, panel_letter, out_name) {
  ds <- datasets[[key]]
  dt <- read_coloc_dir(paste0(ds$dir, ref_suffix))
  dt <- dt[!is.na(PP.H3.abf) & !is.na(PP.H4.abf)]
  cat(ds$label, "(", ref_label, "): ", nrow(dt), "gene x locus tests with computed PP\n")

  p <- ggplot(dt, aes(x = PP.H3.abf, y = PP.H4.abf)) +
    geom_point(alpha = 0.3, size = 1.2, color = ds$color) +
    geom_hline(yintercept = 0.7, linetype = "dashed", color = "firebrick", linewidth = 0.4) +
    geom_hline(yintercept = 0.5, linetype = "dotted", color = "firebrick", linewidth = 0.4) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "PP.H3 (distinct causal variants)",
         y = "PP.H4 (shared causal variant)",
         title = panel_letter) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(size = 14, face = "bold"),
          plot.title.position = "plot")

  ggsave(file.path(base, out_name), p, width = 5, height = 5, dpi = 300, bg = "white")
  cat("Saved:", file.path(base, out_name), "\n\n")
}

panel_map <- list(
  list(key = "primary",    ref_suffix = "",        ref_label = "eQTLGen", panel_letter = "a"),
  list(key = "primary",    ref_suffix = "_gtexlv",  ref_label = "GTEx-LV", panel_letter = "b"),
  list(key = "kramarenko", ref_suffix = "",         ref_label = "eQTLGen", panel_letter = "c"),
  list(key = "kramarenko", ref_suffix = "_gtexlv",  ref_label = "GTEx-LV", panel_letter = "d"),
  list(key = "hcm",        ref_suffix = "",         ref_label = "eQTLGen", panel_letter = "e"),
  list(key = "hcm",        ref_suffix = "_gtexlv",  ref_label = "GTEx-LV", panel_letter = "f"),
  list(key = "dcm",        ref_suffix = "",         ref_label = "eQTLGen", panel_letter = "g"),
  list(key = "dcm",        ref_suffix = "_gtexlv",  ref_label = "GTEx-LV", panel_letter = "h")
)

for (pm in panel_map) {
  ds_label <- tolower(datasets[[pm$key]]$label)
  ref_tag  <- tolower(gsub("-", "", pm$ref_label))
  out_name <- paste0("coloc_quadrant_", ds_label, "_", ref_tag, ".png")
  make_quadrant_plot(pm$key, pm$ref_suffix, pm$ref_label, pm$panel_letter, out_name)
}
