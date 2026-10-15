# 06_CellChat.R

library(Seurat)
library(CellChat)
library(dplyr)
library(tidyr)
library(ggplot2)

seu <- readRDS("data/processed_seurat.rds")

cell_types <- c(
  "Epithelial", "Endothelial", "CD4", "CD8", "Innate lymphoid",
  "B cells", "DC", "Macrophage", "Mono", "Fibroblast", "Mast cells"
)
excluded_states <- c("CD4 Tm", "CD4 Tfh", "CD4 Treg")
focus_pathways <- c("MIF", "CD45", "ADGRE5", "MHC-II", "ITGB2")

meta <- seu[[]] %>%
  tibble::rownames_to_column("cell") %>%
  transmute(
    cell,
    patient_id = as.character(sample_id),
    response = factor(response, levels = c("CR", "non-CR")),
    cell_type = factor(as.character(cell_type), levels = cell_types),
    cell_state = as.character(cell_state)
  ) %>%
  filter(
    !is.na(patient_id), !is.na(response), !is.na(cell_type),
    !is.na(cell_state), !cell_state %in% excluded_states
  )

cell_counts <- meta %>%
  count(response, cell_type, .drop = FALSE, name = "n")

target_counts <- cell_counts %>%
  group_by(cell_type) %>%
  summarise(n_target = min(n), .groups = "drop")

set.seed(2026)

balanced_meta <- meta %>%
  left_join(target_counts, by = "cell_type") %>%
  group_by(response, cell_type) %>%
  group_modify(~ {
    n_take <- unique(.x$n_target)
    .x[sample.int(nrow(.x), size = n_take), , drop = FALSE]
  }) %>%
  ungroup()

balanced_counts <- balanced_meta %>%
  count(response, cell_type, name = "n") %>%
  pivot_wider(names_from = response, values_from = n)

seu_cc <- subset(seu, cells = balanced_meta$cell)
DefaultAssay(seu_cc) <- "RNA"

seu_cc <- NormalizeData(
  seu_cc,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

run_cellchat <- function(x) {
  metadata <- x[[]]
  metadata$cell_type <- factor(
    as.character(metadata$cell_type),
    levels = cell_types
  )

  expression <- GetAssayData(
    x, assay = "RNA", layer = "data"
  )

  cc <- createCellChat(
    object = expression,
    meta = metadata,
    group.by = "cell_type"
  )

  cc@DB <- CellChatDB.human
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, population.size = FALSE)
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)

  cc
}

cellchat_CR <- run_cellchat(
  subset(
    seu_cc,
    cells = balanced_meta$cell[balanced_meta$response == "CR"]
  )
)

cellchat_nonCR <- run_cellchat(
  subset(
    seu_cc,
    cells = balanced_meta$cell[balanced_meta$response == "non-CR"]
  )
)

cellchat_merged <- mergeCellChat(
  list(CR = cellchat_CR, `non-CR` = cellchat_nonCR),
  add.names = c("CR", "non-CR"),
  cell.prefix = TRUE
)

get_pathway_flow <- function(cc) {
  data.frame(
    pathway = cc@netP$pathways,
    information_flow = as.numeric(
      apply(cc@netP$prob, 3, sum, na.rm = TRUE)
    )
  )
}

pathway_flow <- full_join(
  get_pathway_flow(cellchat_CR) %>%
    rename(CR = information_flow),
  get_pathway_flow(cellchat_nonCR) %>%
    rename(`non-CR` = information_flow),
  by = "pathway"
) %>%
  mutate(
    CR = replace_na(CR, 0),
    `non-CR` = replace_na(`non-CR`, 0),
    delta = `non-CR` - CR
  )

focus_flow <- tibble(pathway = focus_pathways) %>%
  left_join(pathway_flow, by = "pathway") %>%
  mutate(
    CR = replace_na(CR, 0),
    `non-CR` = replace_na(`non-CR`, 0),
    delta = `non-CR` - CR
  )

p_count <- compareInteractions(
  cellchat_merged,
  show.legend = FALSE,
  group = c(1, 2)
)

p_weight <- compareInteractions(
  cellchat_merged,
  show.legend = FALSE,
  group = c(1, 2),
  measure = "weight"
)

available_pathways <- intersect(
  focus_pathways,
  union(
    cellchat_CR@netP$pathways,
    cellchat_nonCR@netP$pathways
  )
)

print(focus_flow)
print(p_count)
print(p_weight)

if (length(available_pathways) > 0) {
  p_rank <- rankNet(
    cellchat_merged,
    mode = "comparison",
    signaling = available_pathways,
    color.use = c("#5bacdb", "#ea7fa3"),
    stacked = TRUE,
    do.stat = TRUE,
    comparison = c(1, 2)
  )

  p_outgoing <- netVisual_bubble(
    cellchat_merged,
    sources.use = "CD4",
    signaling = available_pathways,
    comparison = c(1, 2),
    angle.x = 45
  )

  p_incoming <- netVisual_bubble(
    cellchat_merged,
    targets.use = "CD4",
    signaling = available_pathways,
    comparison = c(1, 2),
    angle.x = 45
  )

  print(p_rank)
  print(p_outgoing)
  print(p_incoming)
}
