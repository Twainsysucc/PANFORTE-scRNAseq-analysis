# 02_Seurat_Harmony.R


library(Seurat)
library(harmony)
library(ggplot2)

seu <- readRDS("data/annotated_seurat.rds")
DefaultAssay(seu) <- "RNA"

mt_fraction <- PercentageFeatureSet(seu, pattern = "^MT-") / 100
mt_cutoff <- ifelse(seu$cell_type == "Epithelial", 0.70, 0.20)
seu <- subset(seu, cells = colnames(seu)[!is.na(mt_cutoff) & mt_fraction < mt_cutoff])

seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, nfeatures = 2000, verbose = FALSE)
seu <- ScaleData(seu, features = VariableFeatures(seu), verbose = FALSE)
seu <- RunPCA(seu, features = VariableFeatures(seu), npcs = 30, verbose = FALSE)
seu <- RunHarmony(
  seu,
  group.by.vars = "sample_id",
  reduction.use = "pca",
  reduction.save = "harmony",
  verbose = FALSE
)
seu <- RunUMAP(
  seu,
  reduction = "harmony",
  dims = 1:30,
  reduction.name = "harmony_umap",
  verbose = FALSE
)

umap_df <- as.data.frame(Embeddings(seu, reduction = "harmony_umap")[, 1:2])
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
meta <- seu@meta.data[rownames(umap_df), c("cell_type", "response"), drop = FALSE]
umap_df$cell_type <- meta$cell_type
umap_df$response <- meta$response

cell_types <- sort(unique(as.character(seu$cell_type)))
celltype_colors <- setNames(
  grDevices::hcl.colors(length(cell_types), palette = "Dark 3"),
  cell_types
)

p_umap <- ggplot(umap_df, aes(UMAP_1, UMAP_2)) +
  geom_point(aes(color = cell_type), size = 0.1, alpha = 0.7) +
  scale_color_manual(values = celltype_colors) +
  theme_void(base_size = 16) +
  theme(legend.title = element_blank())

markers <- c(
  "KRT8", "EPCAM", "CD24",
  "PECAM1", "ENG", "VWF",
  "CD3D", "CD4", "CD8A", "CD8B", "KLRC2",
  "CD79A", "MS4A1", "MZB1",
  "CCR7", "CD1C", "CD68", "CD14", "FCGR3A", "SELL",
  "COL3A1", "DCN", "ACTA2",
  "TPSAB1", "CPA3", "KIT"
)

p_markers <- DotPlot(seu, features = markers, group.by = "cell_type") +
  scale_color_gradient2(low = "#5bacdb", mid = "white", high = "#ea7fa3") +
  RotatedAxis() +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

p_density <- ggplot(umap_df, aes(UMAP_1, UMAP_2)) +
  geom_point(aes(color = cell_type), size = 0.1, alpha = 0.5) +
  stat_density_2d(
    aes(fill = after_stat(level), group = response),
    geom = "polygon", bins = 6, alpha = 0.7
  ) +
  stat_density_2d(
    aes(group = response),
    geom = "path", bins = 6, color = "black", linewidth = 0.5
  ) +
  scale_color_manual(name = "Cell Type", values = celltype_colors) +
  scale_fill_gradientn(
    name = "Density",
    colors = c("white", "#ea7fa3", "darkred"),
    limits = c(0, 0.010),
    oob = scales::squish
  ) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_bw() +
  theme(panel.grid = element_blank())

p_composition <- ggplot(
  composition_df,
  aes(x = response, y = proportion, fill = cell_type)
) +
  geom_col(position = "fill", width = 0.6, alpha = 0.8) +
  scale_fill_manual(values = celltype_colors) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "Mean patient-level proportion") +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_umap)
print(p_markers)
print(p_density)
print(p_composition)

