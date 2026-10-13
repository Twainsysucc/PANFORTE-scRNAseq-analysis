# 04_MiloR.R

library(Seurat)
library(miloR)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)

seu <- readRDS("data/processed_seurat.rds")
seu$patient_id <- as.character(seu$sample_id)
seu$response <- factor(seu$response, levels = c("CR", "non-CR"))

run_milo <- function(seu, selected_types) {
  cells <- colnames(seu)[as.character(seu$cell_type) %in% selected_types]
  sub_seu <- subset(seu, cells = cells)
  DefaultAssay(sub_seu) <- "RNA"

  sub_seu <- NormalizeData(sub_seu, verbose = FALSE)
  sub_seu <- FindVariableFeatures(sub_seu, nfeatures = 2000, verbose = FALSE)
  sub_seu <- ScaleData(sub_seu, features = VariableFeatures(sub_seu), verbose = FALSE)
  sub_seu <- RunPCA(
    sub_seu,
    features = VariableFeatures(sub_seu),
    npcs = 30,
    reduction.name = "milo_pca",
    reduction.key = "MILOPC_",
    verbose = FALSE
  )

  metadata <- sub_seu[[]]

  sce <- SingleCellExperiment(
    assays = list(counts = Matrix::Matrix(
      0, nrow = 1, ncol = ncol(sub_seu), sparse = TRUE,
      dimnames = list("dummy", colnames(sub_seu))
    )),
    colData = S4Vectors::DataFrame(metadata)
  )

  reducedDim(sce, "PCA") <- Embeddings(sub_seu, "milo_pca")[colnames(sub_seu), 1:30]
  reducedDim(sce, "UMAP") <- Embeddings(sub_seu, "harmony_umap")[colnames(sub_seu), 1:2]

  milo <- Milo(sce)
  milo <- buildGraph(milo, k = 30, d = 30, reduced.dim = "PCA")

  set.seed(20260913)
  milo <- makeNhoods(
    milo, prop = 0.05, k = 30, d = 30,
    refined = TRUE, reduced_dims = "PCA"
  )

  milo <- countCells(milo, meta.data = metadata, samples = "patient_id")

  design <- unique(metadata[, c("patient_id", "response")])
  stopifnot(!anyDuplicated(design$patient_id))
  rownames(design) <- design$patient_id
  design <- design[colnames(nhoodCounts(milo)), , drop = FALSE]
  stopifnot(!anyNA(design$patient_id))

  milo <- calcNhoodDistance(milo, d = 30, reduced.dim = "PCA")
  da <- testNhoods(
    milo,
    design = ~ response,
    design.df = design,
    norm.method = "logMS",
    fdr.weighting = "graph-overlap",
    robust = TRUE
  )

  da <- annotateNhoods(milo, da, coldata_col = "cell_type")
  da <- annotateNhoods(milo, da, coldata_col = "cell_state")

  da <- da %>%
    mutate(
      state = ifelse(cell_state_fraction >= 0.60, as.character(cell_state), "Mixed"),
      direction = ifelse(logFC > 0, "non-CR enriched", "CR enriched")
    )

  significant <- da %>%
    filter(!is.na(SpatialFDR), SpatialFDR < 0.10)

  summary <- significant %>%
    group_by(state, direction) %>%
    summarise(
      n_neighborhoods = n(),
      median_logFC = median(logFC),
      min_SpatialFDR = min(SpatialFDR),
      .groups = "drop"
    )

  milo <- buildNhoodGraph(milo)

  p_umap <- plotNhoodGraphDA(
    milo, da, alpha = 0.10, layout = "UMAP"
  )

  p_beeswarm <- plotDAbeeswarm(
    filter(da, state != "Mixed"),
    group.by = "state",
    alpha = 0.10
  ) +
    scale_color_gradient2(
      low = "#5bacdb",
      mid = "white",
      high = "#ea7fa3",
      midpoint = 0
    ) +
    theme(panel.grid = element_blank())

  list(
    milo = milo,
    da = da,
    summary = summary,
    p_umap = p_umap,
    p_beeswarm = p_beeswarm
  )
}

epithelial <- run_milo(seu, "Epithelial")
t_cells <- run_milo(seu, c("CD4", "CD8"))

print(epithelial$summary)
print(t_cells$summary)

print(epithelial$p_umap)
print(epithelial$p_beeswarm)
print(t_cells$p_umap)
print(t_cells$p_beeswarm)
