# 03_Cell_composition.R

library(Seurat)
library(dplyr)
library(tidyr)
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

cell_types <- sort(unique(as.character(meta$cell_type)))

composition_df <- meta %>%
  count(patient_id, response, cell_type, name = "n_cells") %>%
  group_by(patient_id, response) %>%
  complete(cell_type = cell_types, fill = list(n_cells = 0L)) %>%
  mutate(proportion = n_cells / sum(n_cells)) %>%
  ungroup()

composition_summary <- composition_df %>%
  group_by(response, cell_type) %>%
  summarise(mean_proportion = mean(proportion), .groups = "drop")

prop_list <- speckle::getTransformedProps(
  clusters = factor(meta$cell_type),
  sample = factor(meta$patient_id),
  transform = "asin"
)
design_df <- as.data.frame(patient_info)
design_df <- design_df[
  match(colnames(prop_list$TransformedProps), design_df$patient_id),
  , drop = FALSE
]
rownames(design_df) <- design_df$patient_id
design_df$group <- factor(
  design_df$response,
  levels = c("CR", "non-CR"),
  labels = c("CR", "nonCR")
)
design_df$batch <- droplevels(factor(design_df$batch))
design_matrix <- model.matrix(
  ~ 0 + group + batch,
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

