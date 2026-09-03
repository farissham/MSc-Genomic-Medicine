library(ggplot2)
library(dplyr)
library(scales)

rg_file <- "/Users/USER-PC/Desktop/hcm_dcm_output/results_ldsc_h2_fix/ldsc/cohort.rg_all_pairs.tsv"
raw <- read.delim(rg_file, stringsAsFactors = FALSE)

id_to_label <- c(
  "hcm_dcm"            = "Primary",
  "kramarenko_ccmtag"  = "Reconstructed",
  "tadros_hcm"         = "HCM",
  "zheng_dcm"          = "DCM"
)

raw$Dataset_X <- id_to_label[raw$id1]
raw$Dataset_Y <- id_to_label[raw$id2]

parse_self_check <- function(log_path) {
  lines <- readLines(log_path)
  hdr_i <- grep("^Summary of Genetic Correlation Results", lines)
  tbl_block <- lines[(hdr_i[1] + 1):length(lines)]
  first_blank <- which(!nzchar(trimws(tbl_block)))[1]
  if (!is.na(first_blank)) tbl_block <- tbl_block[seq_len(first_blank - 1)]
  tbl <- read.table(text = paste(tbl_block, collapse = "\n"), header = TRUE, stringsAsFactors = FALSE)
  list(rg = tbl$rg[1], p = tbl$p[1])
}

self_check_dir <- "/Users/USER-PC/Desktop/hcm_dcm_output/results_comparison_all4/ldsc_self_check"
self_check_ids <- c("hcm_dcm", "kramarenko_ccmtag", "tadros_hcm", "zheng_dcm")

diagonal <- do.call(rbind, lapply(self_check_ids, function(id) {
  res <- parse_self_check(file.path(self_check_dir, paste0(id, "_self_check.log")))
  data.frame(Dataset_X = id_to_label[[id]], Dataset_Y = id_to_label[[id]], rg = res$rg, p = res$p)
}))

mirrored <- data.frame(Dataset_X = raw$Dataset_Y, Dataset_Y = raw$Dataset_X, rg = raw$rg, p = raw$p)

plot_data <- rbind(
  data.frame(Dataset_X = raw$Dataset_X, Dataset_Y = raw$Dataset_Y, rg = raw$rg, p = raw$p),
  mirrored, diagonal
)

x_order <- c("DCM", "HCM", "Reconstructed", "Primary")
y_order_top_to_bottom <- c("Primary", "Reconstructed", "HCM", "DCM")

plot_data$Dataset_X <- factor(plot_data$Dataset_X, levels = x_order)
plot_data$Dataset_Y <- factor(plot_data$Dataset_Y, levels = rev(y_order_top_to_bottom))

plot_data$label <- ifelse(plot_data$p < 0.05,
                          paste0(round(plot_data$rg, 2), "*"),
                          as.character(round(plot_data$rg, 2)))

ggplot(plot_data, aes(x = Dataset_X, y = Dataset_Y, fill = rg)) +
  geom_tile(color = "NA", linewidth = 0.8) +
  geom_text(aes(label = label), color = "black", size = 4.5) +
  scale_fill_gradient2(
    low = "#3498db", mid = "white", high = "#e74c3c",
    midpoint = 0, limits = c(-1, 1), oob = scales::squish,
    name = expression(bold(r[g]))
  ) +
  scale_x_discrete(limits = x_order, drop = FALSE) +
  scale_y_discrete(limits = rev(y_order_top_to_bottom), drop = FALSE) +
  coord_fixed(ratio = 1) +
  labs(x = NULL, y = NULL, title = "a") +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(face = "bold", color = "black"),
    axis.text.y = element_text(face = "bold", color = "black"),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0, size = 14, margin = margin(l = -5)),
    plot.title.position = "plot",
    legend.title = element_text(face = "bold"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave("/Users/USER-PC/Desktop/hcm_dcm_output/ldsc_rg_heatmap.png",
       width = 6.5, height = 6, dpi = 300, bg = "white")