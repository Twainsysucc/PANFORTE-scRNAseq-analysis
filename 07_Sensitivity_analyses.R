# 07_Sensitivity_analyses.R

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(miloR)
library(SingleCellExperiment)
library(CellChat)
library(MatchIt)
library(speckle)
library(limma)

seu <- readRDS("data/processed_seurat.rds")
seu_decontx <- readRDS("data/decontx_corrected_seurat.rds")
clinical_df <- read.csv("data/clinical_metadata.csv", stringsAsFactors = FALSE)
regrowth_ids <- readLines("data/regrowth_patient_ids.txt", warn = FALSE)

DefaultAssay(seu) <- "RNA"
DefaultAssay(seu_decontx) <- "RNA"

seu$patient_id <- as.character(seu$sample_id)
seu$response <- factor(seu$response, levels = c("CR", "non-CR"))

common_cells <- intersect(colnames(seu), colnames(seu_decontx))
stopifnot(setequal(common_cells, colnames(seu)))
seu_decontx <- subset(seu_decontx, cells = common_cells)
for (field in c("patient_id", "response", "cell_type", "cell_state")) {
  seu_decontx@meta.data[[field]] <- seu@meta.data[colnames(seu_decontx), field]
}

mito_fraction <- PercentageFeatureSet(seu, pattern = "^MT-") / 100
cells_mito20 <- colnames(seu)[
  as.character(seu$cell_type) != "Epithelial" | mito_fraction < 0.20
]
cells_ambient30 <- colnames(seu_decontx)[
  !is.na(seu_decontx$decontX_contamination) &
    seu_decontx$decontX_contamination <= 0.30
]

qc_objects <- list(
  Q1 = seu,
  Q2 = subset(seu_decontx, cells = cells_ambient30),
  Q3 = subset(seu, cells = cells_mito20),
  Q4 = subset(seu_decontx, cells = intersect(cells_mito20, cells_ambient30))
)

qc_summary <- tibble(
  QC = names(qc_objects),
  n_cells = vapply(qc_objects, ncol, integer(1))
)

patient_info <- seu[[]] %>%
  distinct(patient_id, response)

regrowth_groups <- list(
  G1 = patient_info,
  G2 = filter(patient_info, !patient_id %in% regrowth_ids),
  G3 = patient_info %>%
    mutate(response = factor(
      if_else(patient_id %in% regrowth_ids, "non-CR", as.character(response)),
      levels = c("CR", "non-CR")
    ))
)

clinical_match <- clinical_df %>%
  select(patient_id, cMRF, cEMVI, T4, N2) %>%
  inner_join(patient_info, by = "patient_id") %>%
  mutate(
    cMRF = factor(cMRF, levels = c(0, 1)),
    cEMVI = factor(cEMVI, levels = c(0, 1)),
    T4 = as.integer(T4),
    N2 = as.integer(N2),
    nonCR = as.integer(response == "non-CR")
  )

set.seed(2026)
match_fit <- MatchIt::matchit(
  nonCR ~ cMRF + cEMVI + T4 + N2,
  data = clinical_match,
  method = "nearest",
  distance = "mahalanobis",
  estimand = "ATT",
  ratio = 1,
  replace = FALSE
)
matched_patients <- MatchIt::match.data(match_fit) %>%
  select(patient_id, response, subclass)
regrowth_groups$G4 <- select(matched_patients, patient_id, response)

sample_summary <- bind_rows(lapply(names(regrowth_groups), function(g) {
  regrowth_groups[[g]] %>%
    count(response, name = "n_patients") %>%
    mutate(group = g)
})) %>% select(group, response, n_patients)


composition_sensitivity <- bind_rows(lapply(c("G1", "G2", "G3"), function(g) {
  labels <- regrowth_groups[[g]]
  dat <- seu[[]] %>%
    select(patient_id, cell_type) %>%
    inner_join(labels, by = "patient_id")
  result <- speckle::propeller(
    clusters = factor(dat$cell_type),
    sample = factor(dat$patient_id),
    group = dat$response,
    transform = "asin",
    robust = FALSE,
    trend = FALSE
  )
  as.data.frame(result) %>%
    tibble::rownames_to_column("cell_type") %>%
    mutate(analysis = g)
}))

matched_meta <- seu[[]] %>%
  select(patient_id, cell_type) %>%
  inner_join(matched_patients, by = "patient_id")
matched_props <- speckle::getTransformedProps(
  clusters = factor(matched_meta$cell_type),
  sample = factor(matched_meta$patient_id),
  transform = "asin"
)
matched_design <- matched_patients[
  match(colnames(matched_props$TransformedProps), matched_patients$patient_id),
]
matched_design$group <- factor(
  matched_design$response,
  levels = c("CR", "non-CR"),
  labels = c("CR", "nonCR")
)
matched_design$subclass <- factor(matched_design$subclass)
paired_matrix <- model.matrix(~ 0 + group + subclass, data = matched_design)
paired_contrast <- limma::makeContrasts(groupnonCR - groupCR, levels = paired_matrix)
composition_G4_paired <- speckle::propeller.ttest(
  prop.list = matched_props,
  design = paired_matrix,
  contrasts = as.numeric(paired_contrast),
  robust = FALSE,
  trend = FALSE,
  sort = TRUE
)

fit_milo <- function(x, selected_types, seed = 20261005) {
  x <- subset(x, cells = colnames(x)[as.character(x$cell_type) %in% selected_types])
  DefaultAssay(x) <- "RNA"
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, nfeatures = 2000, verbose = FALSE)
  x <- ScaleData(x, features = VariableFeatures(x), verbose = FALSE)
  x <- RunPCA(x, npcs = 30, features = VariableFeatures(x), verbose = FALSE)

  metadata <- x[[]]
  sce <- SingleCellExperiment(
    assays = list(counts = Matrix::Matrix(
      0, nrow = 1, ncol = ncol(x), sparse = TRUE,
      dimnames = list("dummy", colnames(x))
    )),
    colData = S4Vectors::DataFrame(metadata)
  )
  reducedDim(sce, "PCA") <- Embeddings(x, "pca")[colnames(x), 1:30, drop = FALSE]

  m <- miloR::Milo(sce)
  m <- miloR::buildGraph(m, k = 30, d = 30, reduced.dim = "PCA")
  set.seed(seed)
  m <- miloR::makeNhoods(
    m, prop = 0.05, k = 30, d = 30,
    refined = TRUE, reduced_dims = "PCA"
  )
  m <- miloR::countCells(m, meta.data = metadata, samples = "patient_id")
  m <- miloR::calcNhoodDistance(m, d = 30, reduced.dim = "PCA")
  m
}

test_milo <- function(m, labels) {
  stopifnot(!anyDuplicated(labels$patient_id))
  counts <- miloR::nhoodCounts(m)
  patients <- intersect(colnames(counts), labels$patient_id)
  miloR::nhoodCounts(m) <- counts[, patients, drop = FALSE]
  design <- labels[match(patients, labels$patient_id), , drop = FALSE]
  design$response <- factor(design$response, levels = c("CR", "non-CR"))
  rownames(design) <- design$patient_id

  da <- miloR::testNhoods(
    m, design = ~ response, design.df = design,
    norm.method = "logMS", fdr.weighting = "graph-overlap", robust = TRUE
  )
  da <- miloR::annotateNhoods(m, da, coldata_col = "cell_state")
  da %>% mutate(
    state = if_else(cell_state_fraction >= 0.60,
                    as.character(cell_state), "Mixed")
  )
}

summarise_epi_milo <- function(da, analysis) {
  tibble(
    analysis = analysis,
    n_neighborhoods = nrow(da),
    n_DA = sum(da$SpatialFDR < 0.10, na.rm = TRUE),
    n_C2_DA = sum(da$state == "C2" & da$SpatialFDR < 0.10, na.rm = TRUE),
    C2_DA_pct = 100 * n_C2_DA / n_neighborhoods
  )
}

milo_epi_qc <- lapply(qc_objects, fit_milo, selected_types = "Epithelial")
epi_qc_summary <- bind_rows(lapply(names(milo_epi_qc), function(q) {
  m <- milo_epi_qc[[q]]
  patients <- qc_objects[[q]][[]] %>% distinct(patient_id, response)
  summarise_epi_milo(test_milo(m, patients), q)
}))

epi_regrowth_milo <- milo_epi_qc$Q1
epi_regrowth_summary <- bind_rows(lapply(names(regrowth_groups), function(g) {
  summarise_epi_milo(test_milo(epi_regrowth_milo, regrowth_groups[[g]]), g)
}))

t_states <- c(
  "CD8 Tn", "CD4 Tn", "CD4 Tfh", "CD4 Tm",
  "CD8 Tem", "CD4 Treg", "CD8 Temra/eff"
)

calculate_roe <- function(x, labels, analysis) {
  meta <- x[[]] %>%
    select(patient_id, cell_type, cell_state) %>%
    inner_join(labels, by = "patient_id") %>%
    filter(cell_type %in% c("CD4", "CD8"), cell_state %in% t_states)
  observed <- table(
    factor(meta$cell_state, levels = t_states),
    factor(meta$response, levels = c("CR", "non-CR"))
  )
  expected <- outer(rowSums(observed), colSums(observed)) / sum(observed)
  values <- observed / expected
  as.data.frame(as.table(values), responseName = "RoE") %>%
    rename(cell_state = Var1, response = Var2) %>%
    mutate(analysis = analysis)
}

roe_qc <- bind_rows(lapply(names(qc_objects), function(q) {
  labels <- qc_objects[[q]][[]] %>% distinct(patient_id, response)
  calculate_roe(qc_objects[[q]], labels, q)
}))
roe_regrowth <- bind_rows(lapply(names(regrowth_groups), function(g) {
  calculate_roe(seu, regrowth_groups[[g]], g)
}))

cell_types <- c(
  "Epithelial", "Endothelial", "CD4", "CD8", "Innate lymphoid",
  "B cells", "DC", "Macrophage", "Mono", "Fibroblast", "Mast cells"
)
excluded_states <- c("CD4 Tm", "CD4 Tfh", "CD4 Treg")
focus_pathways <- c("MIF", "CD45", "ADGRE5", "MHC-II", "ITGB2")

prepare_cc_meta <- function(x, labels) {
  x[[]] %>%
    tibble::rownames_to_column("cell") %>%
    select(cell, patient_id, cell_type, cell_state) %>%
    inner_join(labels, by = "patient_id") %>%
    filter(
      cell_type %in% cell_types,
      !is.na(cell_state),
      !cell_state %in% excluded_states
    )
}

cc_meta_qc <- lapply(qc_objects, function(x) {
  prepare_cc_meta(x, x[[]] %>% distinct(patient_id, response))
})

cc_targets <- bind_rows(lapply(names(cc_meta_qc), function(q) {
  cc_meta_qc[[q]] %>%
    count(response, cell_type, name = "n") %>%
    mutate(QC = q)
})) %>%
  complete(
    QC = names(qc_objects),
    response = factor(c("CR", "non-CR"), levels = c("CR", "non-CR")),
    cell_type = cell_types,
    fill = list(n = 0L)
  ) %>%
  group_by(cell_type) %>%
  summarise(n_target = min(n), .groups = "drop")

run_cellchat_pair <- function(x, meta, target = NULL, seed = 2026) {
  if (is.null(target)) {
    target <- meta %>%
      count(response, cell_type, name = "n") %>%
      complete(response = c("CR", "non-CR"), cell_type = cell_types,
               fill = list(n = 0L)) %>%
      group_by(cell_type) %>%
      summarise(n_target = min(n), .groups = "drop")
  }
  set.seed(seed)
  balanced <- meta %>%
    left_join(target, by = "cell_type") %>%
    group_by(response, cell_type) %>%
    group_modify(~ dplyr::slice_sample(.x, n = unique(.x$n_target))) %>%
    ungroup()

  x <- subset(x, cells = balanced$cell)
  DefaultAssay(x) <- "RNA"
  x <- NormalizeData(x, normalization.method = "LogNormalize",
                     scale.factor = 10000, verbose = FALSE)

  fit_one <- function(response_group) {
    chosen <- balanced$cell[balanced$response == response_group]
    sub <- subset(x, cells = chosen)
    expression <- GetAssayData(sub, assay = "RNA", layer = "data")
    md <- sub[[]]
    md$cell_type <- factor(as.character(md$cell_type), levels = cell_types)
    cc <- CellChat::createCellChat(expression, meta = md, group.by = "cell_type")
    cc@DB <- CellChat::CellChatDB.human
    cc <- CellChat::subsetData(cc)
    cc <- CellChat::identifyOverExpressedGenes(cc)
    cc <- CellChat::identifyOverExpressedInteractions(cc)
    cc <- CellChat::computeCommunProb(cc, population.size = FALSE)
    cc <- CellChat::filterCommunication(cc, min.cells = 10)
    cc <- CellChat::computeCommunProbPathway(cc)
    CellChat::aggregateNet(cc)
  }

  models <- list(CR = fit_one("CR"), `non-CR` = fit_one("non-CR"))
  flows <- lapply(models, function(model) {
    setNames(as.numeric(apply(model@netP$prob, 3, sum, na.rm = TRUE)),
             model@netP$pathways)
  })
  tibble(
    pathway = focus_pathways,
    CR = unname(flows$CR[focus_pathways]),
    `non-CR` = unname(flows$`non-CR`[focus_pathways])
  ) %>%
    mutate(across(c(CR, `non-CR`), ~ replace_na(.x, 0)),
           delta = `non-CR` - CR)
}

cc_qc_summary <- bind_rows(lapply(names(qc_objects), function(q) {
  run_cellchat_pair(qc_objects[[q]], cc_meta_qc[[q]], cc_targets) %>%
    mutate(analysis = q)
}))

cc_matched_summary <- run_cellchat_pair(
  seu,
  prepare_cc_meta(seu, regrowth_groups$G4)
) %>% mutate(analysis = "G4")

p_roe_qc <- ggplot(roe_qc, aes(response, cell_state, fill = RoE)) +
  geom_tile(color = "black") +
  geom_text(aes(label = sprintf("%.2f", RoE)), size = 3) +
  facet_wrap(~ analysis, nrow = 1) +
  scale_fill_gradient2(low = "#5bacdb", mid = "white",
                       high = "#ea7fa3", midpoint = 1) +
  labs(x = NULL, y = NULL, fill = "Ro/e") +
  theme_classic()

p_milo_qc <- ggplot(epi_qc_summary, aes(analysis, C2_DA_pct)) +
  geom_col(fill = "#5bacdb", width = 0.65) +
  labs(x = "QC definition", y = "C2 DA neighborhoods (%)") +
  theme_classic()

p_cellchat_qc <- cc_qc_summary %>%
  pivot_longer(c(CR, `non-CR`), names_to = "response", values_to = "information_flow") %>%
  ggplot(aes(analysis, information_flow, fill = response)) +
  geom_col(position = "dodge") +
  facet_wrap(~ pathway, scales = "free_y") +
  scale_fill_manual(values = c(CR = "#5bacdb", `non-CR` = "#ea7fa3")) +
  labs(x = "QC definition", y = "Information flow", fill = NULL) +
  theme_classic()

print(qc_summary)
print(sample_summary)
print(composition_sensitivity)
print(composition_G4_paired)
print(summary(match_fit, standardize = TRUE)$sum.matched)
print(epi_qc_summary)
print(epi_regrowth_summary)
print(cc_qc_summary)
print(cc_matched_summary)
print(p_roe_qc)
print(p_milo_qc)
print(p_cellchat_qc)
