library(ggplot2)
library(dplyr)
library(scales)

conc_file <- "/Users/USER-PC/Desktop/hcm_dcm_output/results_comparison_all4/concordance/concordance_all_pairs.tsv"
raw <- read.delim(conc_file, stringsAsFactors = FALSE)

id_to_label <- c(
  "hcm_dcm"            = "Primary",
  "kramarenko_ccmtag"  = "Reconstructed",
  "tadros_hcm"         = "HCM",
  "zheng_dcm"          = "DCM"
)

raw$Dataset_X <- id_to_label[raw$id1]
raw$Dataset_Y <- id_to_label[raw$id2]

self_check_dir <- "/Users/USER-PC/Desktop/hcm_dcm_output/results_comparison_all4/concordance_self_check"
self_check_ids <- c("hcm_dcm", "kramarenko_ccmtag", "tadros_hcm", "zheng_dcm")

diagonal <- do.call(rbind, lapply(self_check_ids, function(id) {
  d <- read.delim(file.path(self_check_dir, paste0(id, "_self_check_concordance.tsv")),
                   stringsAsFactors = FALSE)
  data.frame(Dataset_X = id_to_label[[id]], Dataset_Y = id_to_label[[id]],
             conc = d$observed_concordance[1], p = d$pvalue[1])
}))

mirrored <- data.frame(Dataset_X = raw$Dataset_Y, Dataset_Y = raw$Dataset_X,
                        conc = raw$observed_concordance, p = raw$pvalue)

plot_data <- rbind(
  data.frame(Dataset_X = raw$Dataset_X, Dataset_Y = raw$Dataset_Y,
             conc = raw$observed_concordance, p = raw$pvalue),
  mirrored, diagonal
)

x_order <- c("DCM", "HCM", "Reconstructed", "Primary")
y_order_top_to_bottom <- c("Primary", "Reconstructed", "HCM", "DCM")

plot_data$Dataset_X <- factor(plot_data$Dataset_X, levels = x_order)
plot_data$Dataset_Y <- factor(plot_data$Dataset_Y, levels = rev(y_order_top_to_bottom))

plot_data$label <- ifelse(plot_data$p < 0.05,
                           paste0(round(plot_data$conc * 100, 1), "%*"),
                           paste0(round(plot_data$conc * 100, 1), "%"))

ggplot(plot_data, aes(x = Dataset_X, y = Dataset_Y, fill = conc)) +
  geom_tile(color = NA) +
  geom_text(aes(label = label), color = "black", size = 4) +
  scale_fill_gradient(
    low = "white", high = "#2a78d6", limits = c(0, 1),
    labels = scales::percent,
    name = "Concordance"
  ) +
  scale_x_discrete(limits = x_order, drop = FALSE) +
  scale_y_discrete(limits = rev(y_order_top_to_bottom), drop = FALSE) +
  coord_fixed(ratio = 1) +
  labs(x = NULL, y = NULL, title = "b") +
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

ggsave("/Users/USER-PC/Desktop/hcm_dcm_output/concordance_heatmap.png",
       width = 6.5, height = 6, dpi = 300, bg = "white")
