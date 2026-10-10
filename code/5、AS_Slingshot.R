suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(slingshot)
  library(dplyr)
  library(Matrix)
})
input_rds <- '../Fig2_AS_final/03_0_subtype_finalRDS/AS_hdWGCNA_res0.2_Subtype_Final.rds'
out_root <- '../Fig2_AS_final/04_1_Slingshot'
tab_dir <- file.path(out_root, 'tables')
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
reduction_name <- 'umap_AS_all_keepAllRegions'
cluster_col <- 'AS_cluster_res0.2'
subtype_col <- 'AS_subtype'
root_subtype <- 'AS_S1_Homeo'
excluded_cluster <- '5'
excluded_subtype <- 'AS_S5_Cilia'
subtype_levels <- c('AS_S1_Homeo','AS_S2_Demy-Reactive','AS_S3_Myelin-Support','AS_S4_IMM','AS_S5_Cilia','AS_S6_Trans','AS_S7_WM-Support','AS_S8_Adaptive','AS_S9_Stress')
plot_levels <- setdiff(subtype_levels, excluded_subtype)
message('Reading input RDS: ', input_rds)
obj <- readRDS(input_rds)
emb <- as.data.frame(Embeddings(obj, reduction_name))
colnames(emb)[1:2] <- c('AS_UMAP_1', 'AS_UMAP_2')
emb$cell_barcode <- rownames(emb)
meta <- obj@meta.data
meta$cell_barcode <- rownames(meta)
required_cols <- c('cell_barcode', cluster_col, subtype_col)
missing_cols <- setdiff(required_cols, colnames(meta))
if (length(missing_cols) > 0) stop('Missing required metadata columns: ', paste(missing_cols, collapse = ', '))
keep_meta_cols <- intersect(c('cell_barcode', cluster_col, subtype_col, 'Sample_id', 'Study', 'diagnosis', 'donor_id', 'lesion_type', 'matter', 'matter_group'), colnames(meta))
df_all <- emb %>%
  left_join(meta[, keep_meta_cols, drop = FALSE], by = 'cell_barcode') %>%
  mutate(
    AS_cluster_res0.2_chr = as.character(.data[[cluster_col]]),
    AS_subtype_chr = as.character(.data[[subtype_col]]),
    excluded_C5 = AS_cluster_res0.2_chr == excluded_cluster | AS_subtype_chr == excluded_subtype
  )
df <- df_all %>%
  filter(!excluded_C5, !is.na(AS_subtype_chr), is.finite(AS_UMAP_1), is.finite(AS_UMAP_2)) %>%
  mutate(AS_subtype = factor(AS_subtype_chr, levels = plot_levels))
if (!root_subtype %in% df$AS_subtype_chr) stop('Root subtype not found after filtering: ', root_subtype)
filter_summary <- data.frame(
  metric = c('input_cells', 'removed_C5_or_AS_S5_Cilia', 'removed_missing_subtype_or_umap', 'used_cells', 'remaining_subtypes'),
  value = c(nrow(df_all), sum(df_all$excluded_C5, na.rm = TRUE), nrow(df_all) - sum(df_all$excluded_C5, na.rm = TRUE) - nrow(df), nrow(df), paste(levels(droplevels(df$AS_subtype)), collapse = '; '))
)
write.csv(filter_summary, file.path(tab_dir, 'slingshot_input_filter_summary.csv'), row.names = FALSE)
message('Running Slingshot on ', nrow(df), ' cells.')
dummy_counts <- Matrix(0, nrow = 1, ncol = nrow(df), sparse = TRUE)
rownames(dummy_counts) <- 'dummy_feature'
colnames(dummy_counts) <- df$cell_barcode
sce <- SingleCellExperiment(assays = list(counts = dummy_counts), colData = df)
reducedDims(sce)$UMAP <- as.matrix(df[, c('AS_UMAP_1', 'AS_UMAP_2')])
sce <- slingshot(sce, clusterLabels = 'AS_subtype', reducedDim = 'UMAP', start.clus = root_subtype, stretch = 0)
pt <- slingPseudotime(sce)
pt_mean <- rowMeans(pt, na.rm = TRUE)
pt_mean[is.nan(pt_mean)] <- NA_real_
df$slingshot_pseudotime <- pt_mean
curves <- slingCurves(sce)
curve_df <- bind_rows(lapply(seq_along(curves), function(i) {
  s <- curves[[i]]$s
  if (!is.null(curves[[i]]$ord)) s <- s[curves[[i]]$ord, , drop = FALSE]
  as.data.frame(s) %>%
    setNames(c('AS_UMAP_1', 'AS_UMAP_2')) %>%
    mutate(lineage = paste0('Lineage_', i), point_index = row_number())
})) %>%
  arrange(lineage, point_index)
lineage_df <- bind_rows(lapply(seq_along(slingLineages(sce)), function(i) {
  data.frame(lineage = paste0('Lineage_', i), order = seq_along(slingLineages(sce)[[i]]), AS_subtype = slingLineages(sce)[[i]])
}))
cell_out_cols <- intersect(c('cell_barcode', 'AS_UMAP_1', 'AS_UMAP_2', 'AS_cluster_res0.2', 'AS_subtype', 'Sample_id', 'Study', 'diagnosis', 'donor_id', 'lesion_type', 'matter', 'matter_group', 'slingshot_pseudotime'), colnames(df))
write.csv(df[, cell_out_cols, drop = FALSE], file.path(tab_dir, 'slingshot_pseudotime_by_cell.csv'), row.names = FALSE)
write.csv(curve_df, file.path(tab_dir, 'slingshot_curve_coordinates.csv'), row.names = FALSE)
write.csv(curve_df, file.path(tab_dir, 'slingshot_raw_curve_coordinates.csv'), row.names = FALSE)
write.csv(lineage_df, file.path(tab_dir, 'slingshot_lineage_subtype_order.csv'), row.names = FALSE)
