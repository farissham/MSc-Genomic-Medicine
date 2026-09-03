library(ggplot2)
library(scales)

base <- "/Users/USER-PC/Desktop/hcm_dcm_output"

datasets <- list(
  primary    = list(par = file.path(base, "results_hcm_dcm/sbayesrc/hcm_dcm.par"),
                     label = "Primary",    bar_color = "#2a78d6", line_color = "#eb6834"),
  kramarenko = list(par = file.path(base, "results_kramarenko/sbayesrc/kramarenko_ccmtag.par"),
                     label = "Reconstructed", bar_color = "#005F73", line_color = "#CA6702"),
  hcm        = list(par = file.path(base, "results_tadros/sbayesrc/tadros_hcm.par"),
                     label = "HCM",        bar_color = "#4a3aa7", line_color = "#e34948"),
  dcm        = list(par = file.path(base, "results_zheng/sbayesrc/zheng_dcm.par"),
                     label = "DCM",        bar_color = "#e87ba4", line_color = "#008300")
)

max_snp <- 35000

for (key in names(datasets)) {

  ds <- datasets[[key]]
  cat("\n==", ds$label, "==\n")

  par_data <- read.delim(ds$par, header = TRUE, stringsAsFactors = FALSE)
  get_val <- function(item) par_data$Mean[par_data$Item == item]
  get_se  <- function(item) par_data$SD[par_data$Item == item]

  data <- data.frame(
    Component = factor(c("Small\n(C2)", "Medium\n(C3)", "Large\n(C4)", "Very Large\n(C5)"),
                        levels = c("Small\n(C2)", "Medium\n(C3)", "Large\n(C4)", "Very Large\n(C5)")),
    Vg     = c(get_val("Vg2"), get_val("Vg3"), get_val("Vg4"), get_val("Vg5")),
    Vg_SE  = c(get_se("Vg2"),  get_se("Vg3"),  get_se("Vg4"),  get_se("Vg5")),
    NumSnp = c(get_val("NumSnp2"), get_val("NumSnp3"), get_val("NumSnp4"), get_val("NumSnp5"))
  )

  bar_label  <- "SNP Count"
  line_label <- "Proportion of Heritability"

  p <- ggplot(data, aes(x = Component)) +
    geom_col(aes(y = NumSnp, fill = bar_label), width = 0.6, alpha = 0.9) +
    geom_line(aes(y = Vg * max_snp, group = 1, color = line_label), linewidth = 1.2) +
    geom_point(aes(y = Vg * max_snp, color = line_label), size = 3) +
    geom_errorbar(aes(ymin = pmax(0, (Vg - Vg_SE) * max_snp), ymax = pmin(max_snp, (Vg + Vg_SE) * max_snp)),
                  width = 0.15, color = ds$line_color, alpha = 0.6) +
    scale_fill_manual(name = NULL, values = setNames(ds$bar_color, bar_label)) +
    scale_color_manual(name = NULL, values = setNames(ds$line_color, line_label)) +
    scale_y_continuous(
      name = "SNP Count",
      labels = scales::comma,
      limits = c(0, max_snp),
      breaks = seq(0, max_snp, by = 5000),
      expand = c(0, 0),
      sec.axis = sec_axis(~ . / max_snp, name = "Proportion of Heritability",
                          breaks = seq(0, 1, by = 0.20))
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.title.y.left  = element_text(color = "black", face = "bold"),
      axis.text.y.left   = element_text(color = "black"),
      axis.title.y.right = element_text(color = "black", face = "bold"),
      axis.text.y.right  = element_text(color = "black"),
      axis.title.x = element_text(face = "bold", margin = margin(t = 12)),
      legend.position = c(1.00, 1.00),
      legend.justification = c("right", "top"),
      legend.background = element_blank(),
      legend.spacing.y = unit(0, "pt"),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    labs(x = "SBayesRC Mixture Components")

  out_png <- file.path(base, paste0("sbayesrc_architecture_", ds$label, ".png"))
  ggsave(out_png, p, width = 9, height = 5.5, dpi = 300, bg = "white")
  cat("Saved:", out_png, "\n")
}
