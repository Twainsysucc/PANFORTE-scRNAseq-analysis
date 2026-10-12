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
    cell_type = as.character(cell_type)
  ) %>%
  filter(!is.na(patient_id), !is.na(response), !is.na(cell_type))

patient_info <- distinct(meta, patient_id, response)

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

propeller_results <- speckle::propeller(
  clusters = factor(meta$cell_type),
  sample = factor(meta$patient_id),
  group = meta$response,
  transform = "asin",
  robust = FALSE,
  trend = FALSE
) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell_type") %>%
  select(cell_type, P_value = P.Value, FDR)

print(composition_summary)
print(propeller_results)
