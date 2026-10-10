library(Seurat)
library(harmony)
library(dplyr)
Will2<-readRDS("/data/work/Online_data/2025_Will_sn/Cohort2_Validation/Will2_rightgene_names_norm_Atlas.rds")
Martin<-readRDS('/data/work/Online_data/2024Martin_ST_Sn_GSE279183/GSE279183_Normalied_sn_atlas.rds')
Fagiani2025<-readRDS('/data/work/Online_data/2025Absinta_2/GSE277430_norm_snAtlas.rds')
if (anyDuplicated(c(colnames(Will2), colnames(Martin), colnames(Fagiani2025)))) {
  stop("Input cell barcodes must be unique across studies before merging.")
}
Will2$study_new <- "Will2"
Martin$study_new <- "Martin"
Fagiani2025$study_new <- "Fagiani2025"
combined_ms <- merge(x = Will2,
                    y = c(Martin, Fagiani2025),
                    add.cell.ids = NULL)
cols_to_keep <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "sample_id",
                 "type_broad", "sex", "age_at_death", "donor_id",
                  "diagnosis", "lesion_type", "origin", "disease",
                 "age", "pathology", "patient", "annotation_confident",
                 "patient_id", "condition", "ms_class", "celltype", "subtype","study_new")
cols_to_keep <- unique(cols_to_keep)
available_cols <- cols_to_keep[cols_to_keep %in% colnames(combined_ms@meta.data)]
combined_ms@meta.data <- combined_ms@meta.data[, available_cols, drop = FALSE]
combined_ms@meta.data$Celltype <- dplyr::coalesce(
  combined_ms@meta.data$celltype,
  combined_ms@meta.data$annotation_confident,
  combined_ms@meta.data$type_broad
)
combined_ms@meta.data$Disease_Status <- dplyr::coalesce(
  combined_ms@meta.data$diagnosis,
  combined_ms@meta.data$pathology,
  combined_ms@meta.data$condition
)
combined_ms@meta.data$Disease_Status <- recode(combined_ms@meta.data$Disease_Status,
  "MS" = "SPMS",
  "CTR" = "Control"
)
combined_ms@meta.data$Sample_id <- NA_character_
matched_will2 <- intersect(colnames(Will2), colnames(combined_ms))
combined_ms@meta.data[matched_will2, "Sample_id"] <- Will2@meta.data$sample_id[match(matched_will2, colnames(Will2))]
matched_martin <- intersect(colnames(Martin), colnames(combined_ms))
combined_ms@meta.data[matched_martin, "Sample_id"] <- paste0("Ma24_", Martin@meta.data$sample_id[match(matched_martin, colnames(Martin))])
matched_fagiani <- intersect(colnames(Fagiani2025), colnames(combined_ms))
combined_ms@meta.data[matched_fagiani, "Sample_id"] <- paste0("Fa25_", Fagiani2025$orig.ident[match(matched_fagiani, colnames(Fagiani2025))])
combined_ms@meta.data$orig.ident <- combined_ms@meta.data$Sample_id
combined_ms@meta.data$Sample_id <- NULL
combined_ms@meta.data$Patient_ID <- NA_character_
matched_will2 <- intersect(colnames(Will2), colnames(combined_ms))
combined_ms@meta.data[matched_will2, "Patient_ID"] <-
  Will2@meta.data$donor_id[match(matched_will2, colnames(Will2))]
matched_fagiani <- intersect(colnames(Fagiani2025), colnames(combined_ms))
combined_ms@meta.data[matched_fagiani, "Patient_ID"] <-
  paste0("Fa25_", Fagiani2025@meta.data$patient[match(matched_fagiani, colnames(Fagiani2025))])
matched_martin <- intersect(colnames(Martin), colnames(combined_ms))
combined_ms@meta.data[matched_martin, "Patient_ID"] <-
  paste0("Ma24_", Martin@meta.data$patient_id[match(matched_martin, colnames(Martin))])
combined_ms@meta.data$Patient_ID[is.na(combined_ms@meta.data$Patient_ID)] <- "Unassigned"
combined_ms@meta.data <- combined_ms@meta.data %>%
  mutate(lesion_type = ifelse(!is.na(pathology) & pathology != "", pathology, lesion_type))
if (anyNA(combined_ms$orig.ident) || any(combined_ms$orig.ident == "")) {
  stop("Sample IDs are missing after metadata merging.")
}
rm(Will2, Martin, Fagiani2025)
gc()
combined_ms[["percent.mt"]] <- PercentageFeatureSet(combined_ms, pattern = "^MT-")
missing_idx <- which(is.na(combined_ms$nFeature_RNA) | combined_ms$nFeature_RNA == 0)
if (length(missing_idx) > 0) {
  missing_counts <- GetAssayData(combined_ms, assay = "RNA", slot = "counts")[, missing_idx]
  combined_ms@meta.data[missing_idx, "nCount_RNA"] <- colSums(missing_counts)
  combined_ms@meta.data[missing_idx, "nFeature_RNA"] <- colSums(missing_counts > 0)
}
combined_ms_filtered <- subset(
  combined_ms,
  subset = nFeature_RNA > 500 &
           nFeature_RNA < 6000
)
combined_ms_filtered <- FindVariableFeatures(combined_ms_filtered, selection.method = "vst", nfeatures = 2000)
combined_ms_filtered <- ScaleData(combined_ms_filtered)
combined_ms_filtered <- RunPCA(combined_ms_filtered, npcs = 50, verbose = FALSE)
combined_ms_filtered <- RunHarmony(
  combined_ms_filtered,
  group.by.vars = "orig.ident",
  dims.use = 1:50,
  max_iter = 20,
  early_stop = FALSE,
  theta = 1.2,
  sigma = 0.08,
  nclust = 150
)
combined_ms_filtered <- RunUMAP(
  combined_ms_filtered,
  reduction = "harmony",
  dims = 1:50,
  n.neighbors = 20,
  min.dist = 0.1
)
combined_ms_filtered <- FindNeighbors(
  combined_ms_filtered,
  reduction = "harmony",
  dims = 1:25,
  k.param = 20
)
combined_ms_filtered <- FindClusters(
  combined_ms_filtered,
  resolution = 0.7,
  algorithm = 4,
  method = "igraph"
)
combined_ms <- combined_ms_filtered
celltype_groups <- list(
  OL = c("1", "2", "12", "9", "30", "28", "14", "22", "6", "7"),
  AS = c("5", "10", "23", "26", "29"), MG = "4", Vas = c("17", "25"),
  IMM = "24", OPC = "8", NEU = c("21", "3", "16", "11", "13", "27", "18", "20", "19", "15")
)
cluster_to_celltype <- unlist(lapply(names(celltype_groups), function(x) {
  setNames(rep(x, length(celltype_groups[[x]])), celltype_groups[[x]])
}))
combined_ms$celltype_anno <- unname(cluster_to_celltype[as.character(combined_ms$RNA_snn_res.0.7)])
if (anyNA(combined_ms$celltype_anno)) stop("Unmapped Plan C cluster IDs.")
cell_types <- as.character(combined_ms$celltype_anno)
nn_idx <- FNN::get.knn(Embeddings(combined_ms, "umap"), k = 50)$nn.index
offset_status <- rep("Normal", length(cell_types))
for (target_type in unique(cell_types)) {
  target_counts <- rowSums(matrix((cell_types == target_type)[nn_idx], ncol = 50))
  offset_status[cell_types != target_type & target_counts >= 35] <- paste0("Offset_to_", target_type)
}
combined_ms$Offset_Status <- offset_status
rm(nn_idx)
combined_ms$Offset_cell_final_ID <- combined_ms$Offset_Status
if (any(combined_ms$Offset_Status != "Normal")) {
  offset_subset <- subset(combined_ms, subset = Offset_Status != "Normal")
  marker_groups <- list(
    OL = c("PLP1", "MBP", "MOG"), AS = c("GFAP", "SLC1A2", "AQP4"),
    MG = c("CSF1R", "TMEM119", "P2RY12"), Vas = c("CEMIP", "CLDN5", "LAMA2"),
    IMM = c("IGKC", "IGHG1", "FCRL5"), OPC = c("PCDH15", "TNS3", "GPR17"),
    NEU = c("NRGN", "NRG1", "TSHZ2")
  )
  all_features <- unique(unlist(marker_groups))
  if (!all(all_features %in% rownames(offset_subset))) stop("Required annotation markers are missing.")
  offset_subset <- ScaleData(offset_subset, features = all_features)
  exp_matrix <- as.matrix(GetAssayData(offset_subset, slot = "scale.data"))[all_features, , drop = FALSE]
  group_scores <- do.call(cbind, lapply(marker_groups, function(genes) colMeans(exp_matrix[genes, , drop = FALSE])))
  max_group <- colnames(group_scores)[apply(group_scores, 1, which.max)]
  combined_ms@meta.data[colnames(offset_subset), "Offset_cell_final_ID"] <- paste0("Offset_to_", max_group)
  rm(offset_subset, exp_matrix, group_scores)
}
combined_ms$celltype_final <- ifelse(
  combined_ms$Offset_cell_final_ID == "Normal", combined_ms$celltype_anno,
  sub("Offset_to_", "", combined_ms$Offset_cell_final_ID)
)
combined_ms$lesion_type_raw <- combined_ms$lesion_type
combined_ms@meta.data <- combined_ms@meta.data %>% mutate(lesion_type = case_when(
  lesion_type %in% c("CAL", "CA") ~ "CA", lesion_type %in% c("CIL", "CI") ~ "CI",
  lesion_type %in% c("WM", "Ctrl", "control WM") ~ "CtrlWM",
  lesion_type %in% c("A", "AL") ~ "AL", lesion_type == "control cortex" ~ "CtrlGM",
  lesion_type == "demyelinated cortex" ~ "demyGM", lesion_type == "myelinated cortex" ~ "myGM",
  TRUE ~ lesion_type
))
sample_info <- read.csv("../Sample_info.csv", stringsAsFactors = FALSE)
sample_info <- sample_info %>% mutate(across(where(is.character), ~ trimws(., which = "right")))
sample_fields <- c("Study", "lesion_type", "matter", "donor_id", "diagnosis", "sex", "age_at_death")
sample_info <- sample_info %>% select(all_of(c("Sample_ID", sample_fields)))
if (anyNA(sample_info$Sample_ID) || anyDuplicated(sample_info$Sample_ID)) stop("Sample_info must contain unique nonmissing sample IDs.")
sample_idx <- match(combined_ms$orig.ident, sample_info$Sample_ID)
if (anyNA(sample_idx)) stop("Some atlas samples are absent from Sample_info.")
combined_ms@meta.data[sample_fields] <- sample_info[sample_idx, sample_fields, drop = FALSE]
combined_ms@meta.data <- combined_ms@meta.data %>% select(-any_of(c(
  "study_new", "Celltype", "Disease_Status", "Patient_ID",
  "RNA_snn_res.0.4", "RNA_snn_res.0.6", "RNA_snn_res.0.8", "RNA_snn_res.1"
)))
new_colnames <- colnames(combined_ms@meta.data)
new_colnames[new_colnames == "orig.ident"] <- "Sample_id"
new_colnames[new_colnames == "subtype"] <- "subtype_raw"
colnames(combined_ms@meta.data) <- new_colnames
target_order <- c(
  "Sample_id", "Study", "nCount_RNA", "nFeature_RNA", "percent.mt",
  "lesion_type_raw", "subtype_raw", "donor_id", "lesion_type",
  "diagnosis", "matter", "sex", "age_at_death", "RNA_snn_res.0.7",
  "celltype_anno", "Offset_Status", "Offset_cell_final_ID", "celltype_final"
)
combined_ms@meta.data <- combined_ms@meta.data[, c(
  intersect(target_order, colnames(combined_ms@meta.data)), setdiff(colnames(combined_ms@meta.data), target_order)
), drop = FALSE]
out_file <- "../MS_Celltype_meta_final.rds"
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(combined_ms, out_file)
