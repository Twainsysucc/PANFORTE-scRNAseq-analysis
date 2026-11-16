# 05_T_cell_analysis.R

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
