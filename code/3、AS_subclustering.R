suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(tibble)
})
if (!requireNamespace("harmony", quietly = TRUE)) {
  stop("The harmony package is required for this run but is not installed.")
}
base_dir <- "../results_fig3_AS"
as_rds <- "../MS_AS_NicheRisk_Project/data/CombinedMS_AS_allregion.rds"
outdir <- file.path(base_dir, "01_AS_subclustering_GM_SubcorticalWM_WM_aware_final")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(as_rds)) {
  stop("AS subset RDS not found: ", as_rds)
}
as_obj <- readRDS(as_rds)
if (!inherits(as_obj, "Seurat")) {
  stop("Input is not a Seurat object.")
}
cat("Loaded AS object:\n")
print(as_obj)
cat("\n")
required_cols <- c(
  "celltype_final",
  "matter",
  "lesion_type",
  "Sample_id",
  "donor_id",
  "Study",
  "diagnosis"
)
missing_cols <- setdiff(required_cols, colnames(as_obj@meta.data))
if (length(missing_cols) > 0) {
  stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))
}
analysis_name <- "AS_all_keepAllRegions"
analysis_dir <- file.path(outdir, "01_AS_all_keep_GM_SubcorticalWM_WM_Mix")
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
DefaultAssay(as_obj) <- "RNA"
cat("\n============================================================\n")
cat("Running main analysis:", analysis_name, "\n")
cat("Cells:", ncol(as_obj), "\n")
cat("============================================================\n\n")
cat("NormalizeData...\n")
as_obj <- NormalizeData(
  object = as_obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = TRUE
)
cat("FindVariableFeatures...\n")
as_obj <- FindVariableFeatures(
  object = as_obj,
  selection.method = "vst",
  nfeatures = 3000,
  verbose = TRUE
)
variable_genes <- VariableFeatures(as_obj)
write_csv(
  data.frame(variable_gene = variable_genes),
  file.path(analysis_dir, paste0(analysis_name, "_variable_genes_top3000.csv"))
)
cat("ScaleData without vars.to.regress...\n")
as_obj <- ScaleData(
  object = as_obj,
  features = variable_genes,
  verbose = TRUE
)
cat("RunPCA...\n")
as_obj <- RunPCA(
  object = as_obj,
  features = variable_genes,
  npcs = 50,
  verbose = TRUE
)
reduction_use <- "pca"
dims_use <- 1:30
cat("Running Harmony by Study...\n")
as_obj <- harmony::RunHarmony(
  object = as_obj,
  group.by.vars = "Study",
  reduction.use = "pca",
  dims.use = dims_use,
  verbose = TRUE
)
reduction_use <- "harmony"
cat("Using reduction:", reduction_use, "\n")
umap_name <- paste0("umap_", analysis_name)
cat("RunUMAP...\n")
as_obj <- RunUMAP(
  object = as_obj,
  reduction = reduction_use,
  dims = dims_use,
  reduction.name = umap_name,
  reduction.key = paste0("UMAP_", analysis_name, "_"),
  verbose = TRUE
)
cat("FindNeighbors...\n")
as_obj <- FindNeighbors(
  object = as_obj,
  reduction = reduction_use,
  dims = dims_use,
  verbose = TRUE
)
resolutions <- c(0.1, 0.2, 0.3, 0.4, 0.5)
for (res in resolutions) {
  cat("FindClusters resolution:", res, "\n")
  as_obj <- FindClusters(
    object = as_obj,
    resolution = res,
    algorithm = 1,
    verbose = TRUE
  )
  cluster_col_candidates <- c(
    paste0("RNA_snn_res.", res),
    paste0("SCT_snn_res.", res),
    paste0("integrated_snn_res.", res),
    "seurat_clusters"
  )
  cluster_col <- cluster_col_candidates[
    cluster_col_candidates %in% colnames(as_obj@meta.data)
  ][1]
  if (is.na(cluster_col)) {
    stop("Could not find cluster column for resolution ", res)
  }
  final_res_col <- paste0(analysis_name, "_res_", res)
  as_obj@meta.data[[final_res_col]] <- as.character(as_obj@meta.data[[cluster_col]])
  write_csv(
    as_obj@meta.data %>%
      rownames_to_column("cell_barcode") %>%
      select(cell_barcode, all_of(final_res_col)),
    file.path(analysis_dir, paste0(analysis_name, "_cluster_assignment_res", res, ".csv"))
  )
}
saveRDS(
  as_obj,
  file.path(analysis_dir, paste0(analysis_name, "_subclustered.rds"))
)
