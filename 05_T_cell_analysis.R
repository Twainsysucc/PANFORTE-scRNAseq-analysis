# 05_T_cell_analysis.R

# ------------------------------------------------------------
# T cells Ro/e plots
# ------------------------------------------------------------

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(speckle)
library(ComplexHeatmap)
library(circlize)

seu <- readRDS("data/processed_seurat.rds")
meta <- seu[[]] %>%
  transmute(
    patient_id = as.character(sample_id),
    response = factor(response, levels = c("CR", "non-CR")),
    cell_type = as.character(cell_type),
    cell_state = as.character(cell_state)
  ) %>%
  filter(!is.na(patient_id), !is.na(response))
state_map <- tibble(
  cell_state = c(
    "CD8 Tn", "CD4 Tn", "CD4 Tfh", "CD4 Tm",
    "CD8 Tem", "CD4 Treg", "CD8 Temra/eff"
  ),
  cell_type = c("CD8", "CD4", "CD4", "CD4", "CD8", "CD4", "CD8")
)
patient_info <- distinct(meta, patient_id, response)
t_meta <- meta %>%
  filter(cell_type %in% c("CD4", "CD8")) %>%
  inner_join(state_map, by = c("cell_state", "cell_type"))

observed <- table(
  factor(t_meta$cell_state, levels = state_map$cell_state),
  factor(t_meta$response, levels = c("CR", "non-CR"))
)
expected <- outer(rowSums(observed), colSums(observed)) / sum(observed)
roe <- observed / expected
chi_sq <- suppressWarnings(chisq.test(observed, correct = FALSE))
p_values <- 2 * pnorm(-abs(chi_sq$stdres))
fdr_roe <- matrix(
  p.adjust(as.vector(p_values), method = "BH"),
  nrow = nrow(roe),
  dimnames = dimnames(roe)
)
sig_roe <- ifelse(
  fdr_roe < 0.001, "***",
  ifelse(fdr_roe < 0.01, "**", ifelse(fdr_roe < 0.05, "*", ""))
)
p_roe <- ComplexHeatmap::Heatmap(
  as.matrix(roe),
  name = "Ro/e Index",
  col = circlize::colorRamp2(
    c(min(roe), 1, max(roe)),
    c("#5bacdb", "white", "#ea7fa3")
  ),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_names_side = "right",
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid::grid.text(sprintf("%.2f", roe[i, j]), x, y + grid::unit(1.5, "mm"))
    grid::grid.text(sig_roe[i, j], x, y - grid::unit(2, "mm"))
  }
)
print(p_roe)



# ------------------------------------------------------------
# Monocle 2 trajectory plots
# ------------------------------------------------------------

library(monocle)
library(ggplot2)

base_trajectory <- monocle::plot_cell_trajectory(
  trajectory_cds,
  show_branch_points = TRUE,
  color_by = "cell_state",
  cell_size = 0
)

plot_trajectory_group <- function(response_group) {
  group_data <- trajectory_df[
    !is.na(trajectory_df$response) &
      trajectory_df$response == response_group,
    ,
    drop = FALSE
  ]

  base_trajectory +
    geom_point(
      data = group_data,
      aes(x = x, y = y, color = cell_state),
      size = 2.5,
      alpha = 0.9,
      inherit.aes = FALSE
    ) +
    scale_color_manual(
      values = cell_state_colors,
      name = "Cell state",
      drop = FALSE
    ) +
    labs(title = paste("Trajectory –", response_group)) +
    theme(
      panel.border = element_rect(
        color = "black", fill = NA, linewidth = 1
      ),
      plot.background = element_rect(
        color = "black", fill = NA, linewidth = 1
      )
    )
}

p_trajectory_cr <- plot_trajectory_group("CR")
p_trajectory_noncr <- plot_trajectory_group("non-CR")

print(p_trajectory_cr)
print(p_trajectory_noncr)



# ------------------------------------------------------------
# CellChat: outgoing signaling information flow
# ------------------------------------------------------------

library(CellChat)

response_colors <- c(
  "CR" = "#2d83aa",
  "non-CR" = "#bc3a45"
)

pathways_to_plot <- c(
  "MIF", "CD45", "ADGRE5", "MHC-II", "ITGB2"
)

source_cell <- "CD4 Tn"

p_signaling_flow <- CellChat::rankNet(
  cellchat_comparison,
  mode = "comparison",
  sources.use = source_cell,
  signaling = pathways_to_plot,
  color.use = unname(response_colors),
  stacked = TRUE,
  do.stat = TRUE,
  comparison = c(1, 2)
)

print(p_signaling_flow)



# ------------------------------------------------------------
# CellChat: outgoing network from CD4 Tn cells
# ------------------------------------------------------------

source_cell <- "CD4 Tn"
network_title <- "CD4 Tn outgoing signaling"

interaction_matrix <- cellchat_network@net$weight
interaction_matrix <- interaction_matrix[
  cell_order,
  cell_order,
  drop = FALSE
]

outgoing_matrix <- matrix(
  0,
  nrow = nrow(interaction_matrix),
  ncol = ncol(interaction_matrix),
  dimnames = dimnames(interaction_matrix)
)

outgoing_matrix[source_cell, ] <-
  interaction_matrix[source_cell, ]

cell_numbers <- table(cellchat_network@idents)

node_sizes <- as.numeric(
  cell_numbers[rownames(interaction_matrix)]
)
node_colors <- unname(
  cell_type_colors[rownames(interaction_matrix)]
)

par(
  mfrow = c(1, 1),
  xpd = TRUE,
  mar = c(1, 1, 1, 1)
)

CellChat::netVisual_circle(
  outgoing_matrix,
  vertex.weight = node_sizes,
  weight.scale = TRUE,
  arrow.size = 0.05,
  arrow.width = 1,
  edge.weight.max = max(interaction_matrix),
  title.name = network_title,
  color.use = node_colors
)
