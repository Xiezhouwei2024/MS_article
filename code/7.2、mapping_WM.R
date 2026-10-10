suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(spacexr)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(grid)
  library(svglite)
  library(ragg)
  library(yaml)
})
reference_rds <- Sys.getenv(
  "RCTD_REFERENCE_RDS",
  unset = "../MS_AS_NicheRisk_Project/data/RCTD_RefAtlas_Downsample/RCTD_reference_noNEU_2000perType_bySample_exclude_AS_Ma24_MS549H_MS549T.rds"
)
spatial_rds <- Sys.getenv(
  "RCTD_SPATIAL_RDS",
  unset = "/Users/xzwi/Desktop/MS🧠项目/亚型赋予meta的RDS文件/Dr2_subtype_769BR_D6.rds"
)
output_dir <- Sys.getenv(
  "RCTD_OUTPUT_DIR",
  unset = "../Fig2_Stereo_RCTD/Donor2_noNEU_refatlas_RCTD"
)
sample_label <- Sys.getenv("RCTD_SAMPLE_LABEL", unset = "Donor2")
max_cores <- as.integer(Sys.getenv("RCTD_MAX_CORES", unset = "8"))
doublet_mode <- Sys.getenv("RCTD_DOUBLET_MODE", unset = "doublet")
force_run <- tolower(Sys.getenv("RCTD_FORCE", unset = "false")) %in% c("1", "true", "yes")
exclude_region_col <- Sys.getenv("RCTD_EXCLUDE_REGION_COL", unset = "region")
exclude_region_value <- Sys.getenv("RCTD_EXCLUDE_REGION_VALUE", unset = "")
patch_detect_cores_if_needed <- function(fallback_cores) {
  detected <- parallel::detectCores()
  if (!is.na(detected)) {
    return(invisible(detected))
  }
  message("parallel::detectCores() returned NA; using configured fallback max_cores inside this R session.")
  fallback_cores <- as.integer(fallback_cores)
  parallel_namespace <- asNamespace("parallel")
  unlockBinding("detectCores", parallel_namespace)
  assign(
    "detectCores",
    function(logical = TRUE) fallback_cores,
    envir = parallel_namespace
  )
  lockBinding("detectCores", parallel_namespace)
  invisible(fallback_cores)
}
detected_or_fallback_cores <- patch_detect_cores_if_needed(max_cores)
if (is.na(parallel::detectCores())) {
  max_cores <- as.integer(max_cores)
} else if (detected_or_fallback_cores == 1L && max_cores > 1L) {
  max_cores <- 1L
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "objects"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "logs"), recursive = TRUE, showWarnings = FALSE)
rctd_rds <- file.path(output_dir, "objects", paste0(sample_label, "_noNEU_refatlas_celltype_final_RCTD_result.rds"))
annotated_spatial_rds <- file.path(output_dir, "objects", paste0(sample_label, "_spatial_with_RCTD_noNEU_refatlas_celltype_final.rds"))
weights_csv <- file.path(output_dir, "tables", paste0(sample_label, "_RCTD_celltype_final_weights.csv"))
comparison_csv <- file.path(output_dir, "tables", paste0(sample_label, "_RCTD_original_vs_predicted_comparison.csv"))
confusion_csv <- file.path(output_dir, "tables", paste0(sample_label, "_RCTD_compact_confusion_counts.csv"))
confusion_prop_csv <- file.path(output_dir, "tables", paste0(sample_label, "_RCTD_compact_confusion_row_proportions.csv"))
summary_csv <- file.path(output_dir, "tables", paste0(sample_label, "_RCTD_compact_consistency_summary.csv"))
gene_audit_csv <- file.path(output_dir, "tables", paste0(sample_label, "_RCTD_gene_overlap_audit.csv"))
run_summary_md <- file.path(output_dir, paste0(sample_label, "_RCTD_celltype_final_summary.md"))
params_file <- file.path(output_dir, "logs", paste0(sample_label, "_RCTD_celltype_final_plot_parameters.yaml"))
session_file <- file.path(output_dir, "logs", paste0(sample_label, "_RCTD_celltype_final_session_info.txt"))
stopifnot(file.exists(reference_rds), file.exists(spatial_rds))
get_counts <- function(object, assay = DefaultAssay(object)) {
  tryCatch(
    GetAssayData(object, assay = assay, layer = "counts"),
    error = function(e) GetAssayData(object, assay = assay, slot = "counts")
  )
}
count_integer_audit <- function(counts, label) {
  nonzero_values <- counts@x
  distance_to_integer <- abs(nonzero_values - round(nonzero_values))
  data.frame(
    matrix = label,
    nnz_before_round = length(nonzero_values),
    noninteger_nnz_before_round = sum(distance_to_integer > 1e-6),
    noninteger_fraction_before_round = ifelse(length(nonzero_values) > 0, sum(distance_to_integer > 1e-6) / length(nonzero_values), 0),
    max_abs_distance_to_integer = ifelse(length(nonzero_values) > 0, max(distance_to_integer), 0),
    stringsAsFactors = FALSE
  )
}
round_for_rctd <- function(counts, label) {
  audit <- count_integer_audit(counts, label)
  if (audit$noninteger_nnz_before_round > 0) {
    message(
      label,
      " has ",
      audit$noninteger_nnz_before_round,
      " non-integer nonzero values; rounding temporary RCTD input."
    )
    counts@x <- round(counts@x)
    counts <- drop0(counts)
  }
  list(counts = counts, audit = audit)
}
message("Loading reference: ", reference_rds)
reference_object <- readRDS(reference_rds)
if (is.list(reference_object) && "seurat" %in% names(reference_object)) {
  reference_object <- reference_object$seurat
}
stopifnot(inherits(reference_object, "Seurat"))
stopifnot("celltype_final" %in% colnames(reference_object@meta.data))
message("Loading spatial object: ", spatial_rds)
spatial <- readRDS(spatial_rds)
stopifnot(inherits(spatial, "Seurat"))
stopifnot(all(c("coord_x", "coord_y", "celltype") %in% colnames(spatial@meta.data)))
initial_spatial_locations <- ncol(spatial)
excluded_region_locations <- 0L
if (!identical(exclude_region_value, "")) {
  if (!exclude_region_col %in% colnames(spatial@meta.data)) {
    stop("Requested region exclusion column is missing in spatial metadata: ", exclude_region_col)
  }
  keep_region <- is.na(spatial@meta.data[[exclude_region_col]]) |
    spatial@meta.data[[exclude_region_col]] != exclude_region_value
  excluded_region_locations <- sum(!keep_region)
  message(
    "Excluding ",
    excluded_region_locations,
    " spatial locations where ",
    exclude_region_col,
    " == ",
    exclude_region_value
  )
  spatial <- subset(spatial, cells = rownames(spatial@meta.data)[keep_region])
}
ref_counts <- get_counts(reference_object, assay = "RNA")
spatial_counts <- get_counts(spatial, assay = DefaultAssay(spatial))
common_genes <- intersect(rownames(ref_counts), rownames(spatial_counts))
if (length(common_genes) < 500L) {
  stop("Too few common genes for RCTD: ", length(common_genes))
}
gene_audit <- data.frame(
  reference_genes = nrow(ref_counts),
  spatial_genes = nrow(spatial_counts),
  common_genes = length(common_genes),
  reference_cells = ncol(ref_counts),
  spatial_locations = ncol(spatial_counts),
  stringsAsFactors = FALSE
)
write.csv(gene_audit, gene_audit_csv, row.names = FALSE)
ref_counts <- ref_counts[common_genes, , drop = FALSE]
spatial_counts <- spatial_counts[common_genes, , drop = FALSE]
ref_round <- round_for_rctd(ref_counts, "reference")
spatial_round <- round_for_rctd(spatial_counts, "spatial")
ref_counts <- ref_round$counts
spatial_counts <- spatial_round$counts
count_rounding_audit <- bind_rows(ref_round$audit, spatial_round$audit)
gene_audit <- bind_cols(
  gene_audit[rep(1, nrow(count_rounding_audit)), , drop = FALSE],
  count_rounding_audit
)
write.csv(gene_audit, gene_audit_csv, row.names = FALSE)
ref_labels <- factor(reference_object$celltype_final)
names(ref_labels) <- colnames(reference_object)
keep_ref <- !is.na(ref_labels) & colSums(ref_counts) > 0
ref_counts <- ref_counts[, keep_ref, drop = FALSE]
ref_labels <- droplevels(ref_labels[keep_ref])
ref_n_umi <- Matrix::colSums(ref_counts)
spatial_meta <- spatial@meta.data[colnames(spatial_counts), , drop = FALSE]
coords <- spatial_meta[, c("coord_x", "coord_y")]
coords$coord_x <- as.numeric(coords$coord_x)
coords$coord_y <- as.numeric(coords$coord_y)
rownames(coords) <- rownames(spatial_meta)
spatial_n_umi <- Matrix::colSums(spatial_counts)
keep_spatial <- spatial_n_umi > 0 & stats::complete.cases(coords)
spatial_counts <- spatial_counts[, keep_spatial, drop = FALSE]
coords <- coords[keep_spatial, , drop = FALSE]
spatial_meta <- spatial_meta[keep_spatial, , drop = FALSE]
spatial_n_umi <- spatial_n_umi[keep_spatial]
message("Building RCTD reference and puck")
reference <- Reference(
  cell_types = ref_labels,
  counts = ref_counts,
  nUMI = ref_n_umi,
  n_max_cells = ncol(ref_counts)
)
puck <- SpatialRNA(coords, spatial_counts, spatial_n_umi)
if (!file.exists(rctd_rds) || force_run) {
  message("Running RCTD with max_cores=", max_cores, ", doublet_mode=", doublet_mode)
  my_rctd <- create.RCTD(
    spatialRNA = puck,
    reference = reference,
    max_cores = max_cores,
    CELL_MIN_INSTANCE = 5,
    UMI_min = 1,
    UMI_min_sigma = 1
  )
  my_rctd <- run.RCTD(my_rctd, doublet_mode = doublet_mode)
  saveRDS(my_rctd, rctd_rds, compress = "gzip")
} else {
  message("Using existing RCTD result: ", rctd_rds)
  my_rctd <- readRDS(rctd_rds)
}
weights <- as.matrix(my_rctd@results$weights)
row_sums <- rowSums(weights)
norm_weights <- weights
valid_rows <- row_sums > 0
norm_weights[valid_rows, ] <- sweep(weights[valid_rows, , drop = FALSE], 1, row_sums[valid_rows], "/")
norm_weights[!valid_rows, ] <- NA_real_
predicted <- rep(NA_character_, nrow(norm_weights))
predicted[valid_rows] <- colnames(norm_weights)[max.col(norm_weights[valid_rows, , drop = FALSE], ties.method = "first")]
max_weight <- rep(NA_real_, nrow(norm_weights))
max_weight[valid_rows] <- apply(norm_weights[valid_rows, , drop = FALSE], 1, max, na.rm = TRUE)
results_df <- my_rctd@results$results_df
result_extra <- data.frame(barcode = rownames(norm_weights), stringsAsFactors = FALSE)
if (!is.null(results_df) && nrow(results_df) > 0L) {
  tmp <- as.data.frame(results_df)
  tmp$barcode <- rownames(tmp)
  result_extra <- left_join(result_extra, tmp, by = "barcode")
}
weights_out <- data.frame(
  barcode = rownames(norm_weights),
  predicted_celltype_final = predicted,
  max_norm_weight = max_weight,
  norm_weights,
  check.names = FALSE
)
weights_out <- left_join(weights_out, result_extra, by = "barcode")
write.csv(weights_out, weights_csv, row.names = FALSE)
spatial[["RCTD_celltype_final"]] <- NA_character_
spatial[["RCTD_max_norm_weight"]] <- NA_real_
spatial@meta.data[rownames(norm_weights), "RCTD_celltype_final"] <- predicted
spatial@meta.data[rownames(norm_weights), "RCTD_max_norm_weight"] <- max_weight
for (cell_type in colnames(norm_weights)) {
  column_name <- paste0("RCTD_weight_", cell_type)
  spatial[[column_name]] <- NA_real_
  spatial@meta.data[rownames(norm_weights), column_name] <- norm_weights[, cell_type]
}
saveRDS(spatial, annotated_spatial_rds, compress = "gzip")
