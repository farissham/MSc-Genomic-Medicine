library(ggplot2)
library(dplyr)
library(scales)

rg_file <- "/Users/USER-PC/Desktop/hcm_dcm_output/results_comparison_all4_lv/cohort.rg_all_pairs.tsv"
raw <- read.delim(rg_file, stringsAsFactors = FALSE)

disease_ids <- c(
  "hcm_dcm"           = "Primary"
)

trait_ids <- c(
  "lv_lvedvi"     = "LVedv",
  "lv_lvesvi"     = "LVesv",
  "lv_lvmi"       = "LVm",
  "lv_lvef"       = "LVef",
  "lv_lvconc"     = "LVconc",
  "lv_maxwt"      = "MaxWT",
  "lv_meanwt"     = "MeanWT",
  "lv_strainrad"  = "strain-rad",
  "lv_strainlong" = "strain-long",
  "lv_straincirc" = "strain-circ"
)

disease_order <- c("Primary")
trait_order   <- unname(trait_ids)

plot_data <- do.call(rbind, lapply(names(disease_ids), function(d_id) {
  do.call(rbind, lapply(names(trait_ids), function(t_id) {
    hit <- raw[(raw$id1 == d_id & raw$id2 == t_id) | (raw$id1 == t_id & raw$id2 == d_id), ]
    if (nrow(hit) == 0) {
      warning(paste("No rg pair found for", d_id, "vs", t_id))
      return(data.frame(Dataset = disease_ids[[d_id]], Trait = trait_ids[[t_id]], rg = NA, p = NA))
    }
    data.frame(Dataset = disease_ids[[d_id]], Trait = trait_ids[[t_id]], rg = hit$rg[1], p = hit$p[1])
  }))
}))

plot_data$Dataset <- factor(plot_data$Dataset, levels = rev(disease_order))
plot_data$Trait <- factor(plot_data$Trait, levels = trait_order)

plot_data$label <- ifelse(plot_data$p < 0.05,
                           paste0(round(plot_data$rg, 2), "*"),
                           as.character(round(plot_data$rg, 2)))

ggplot(plot_data, aes(x = Trait, y = Dataset, fill = rg)) +
  geom_tile(color = "NA", linewidth = 0.8) +
  geom_text(aes(label = label), color = "black", size = 4) +
  scale_fill_gradient2(
    low = "#3498db", mid = "white", high = "#e74c3c",
    midpoint = 0, limits = c(-1, 1), oob = scales::squish,
    name = expression(bold(r[g]))
  ) +
  scale_x_discrete(limits = trait_order, drop = FALSE) +
  scale_y_discrete(limits = rev(disease_order), drop = FALSE) +
  labs(x = NULL, y = NULL, title = "c") +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(face = "bold", color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(face = "bold", color = "black"),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0, size = 14, margin = margin(l = -5)),
    plot.title.position = "plot",
    legend.title = element_text(face = "bold"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave("/Users/USER-PC/Desktop/hcm_dcm_output/ldsc_rg_lv_heatmap.png",
       width = 10, height = 2.5, dpi = 300, bg = "white")
