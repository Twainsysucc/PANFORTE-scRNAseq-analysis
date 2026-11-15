# 03_Cell_composition.R

library(Seurat)
library(dplyr)
library(speckle)

seu <- readRDS("data/processed_seurat.rds")

meta <- seu[[]] %>%
  transmute(
    patient_id = as.character(sample_id),
    response = factor(response, levels = c("CR", "non-CR")),
    batch = factor(batch),
    cell_type = as.character(cell_type)
  )

patient_info <- distinct(meta, patient_id, response, batch)

prop_list <- speckle::getTransformedProps(
  clusters = factor(meta$cell_type),
  sample = factor(meta$patient_id),
  transform = "asin"
)

sample_order <- colnames(prop_list$TransformedProps)
idx <- match(sample_order, patient_info$patient_id)
design_df <- as.data.frame(patient_info[idx, ])
rownames(design_df) <- design_df$patient_id

design_df$group <- factor(
  design_df$response,
  levels = c("CR", "non-CR"),
  labels = c("CR", "nonCR")
)
design_df$batch <- droplevels(factor(design_df$batch))

design_formula <- if (nlevels(design_df$batch) > 1) {
  ~ 0 + group + batch
} else {
  message("Only one sequencing batch: fitting a group-only model.")
  ~ 0 + group
}

design_matrix <- model.matrix(
  design_formula,
  data = design_df
)

contrast <- limma::makeContrasts(
  groupnonCR - groupCR,
  levels = design_matrix
)

propeller_results <- speckle::propeller.ttest(
  prop.list = prop_list,
  design = design_matrix,
  contrasts = as.numeric(contrast),
  robust = FALSE,
  trend = FALSE,
  sort = TRUE
) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell_type") %>%
  select(cell_type, P_value = P.Value, FDR)

print(propeller_results)
