suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})
if (!exists("combined_ms")) {
  combined_ms <- readRDS("/data/work/Online_data/MS_Celltype_meta_final.rds")
}
if (!"RNA_snn_res.0.7" %in% colnames(combined_ms@meta.data)) {
  stop("RNA_snn_res.0.7 is missing.")
}
output_dir <- "/data/work/Online_data/1_fig1_reproducible_outputs/FigS_cluster_marker_dotplot"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_base <- file.path(output_dir, "FigS_cluster_marker_dotplot")
marker_csv_candidates <- c(
  "/data/work/Online_data/AllMarkers_Res0.7_Sampling3000.csv",
  "AllMarkers_Res0.7_Sampling3000.csv"
)
marker_csv <- marker_csv_candidates[file.exists(marker_csv_candidates)][1]
if (!is.na(marker_csv)) {
  all_markers_fast <- read.csv(marker_csv, row.names = 1)
} else {
  Idents(combined_ms) <- "RNA_snn_res.0.7"
  all_markers_fast <- FindAllMarkers(
    combined_ms,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25,
    test.use = "wilcox",
    max.cells.per.ident = 3000,
    min.diff.pct = 0.1,
    verbose = TRUE
  )
  marker_csv <- "/data/work/Online_data/AllMarkers_Res0.7_Sampling3000.csv"
  write.csv(all_markers_fast, marker_csv)
}
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
})
if (!exists("combined_ms")) {
  combined_ms <- readRDS("/data/work/Online_data/MS_Celltype_meta_final.rds")
}
output_dir <- "/data/work/Online_data/1_fig1_reproducible_outputs/FigS_reference_correlation"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_base <- file.path(output_dir, "FigS_reference_correlation_heatmap")
source_csv_candidates <- c(
  "/data/work/Online_data/Table_S1_CellType_Correlation_Full.csv",
  "Table_S1_CellType_Correlation_Full.csv"
)
source_csv <- source_csv_candidates[file.exists(source_csv_candidates)][1]
if (!is.na(source_csv)) {
  plot_data <- read.csv(source_csv, check.names = FALSE) %>%
    transmute(
      My_Data = My_Atlas_CellType,
      Ref_Data = Reference_Atlas_CellType,
      Correlation = Pearson_Correlation,
      Significance = Significance_Mark
    )
} else {
  if (!exists("ref_seurat")) {
    ref_seurat <- readRDS("/data/work/Online_data/2025_Will_sn/Cohort1/Will_1_right_genenames_norm_Atlas.rds")
  }
  Idents(combined_ms) <- "celltype_final"
  Idents(ref_seurat) <- "type_broad"
  common_genes <- intersect(rownames(combined_ms), rownames(ref_seurat))
  pb_my <- AggregateExpression(
    combined_ms,
    assays = "RNA",
    slot = "counts",
    features = common_genes,
    return.seurat = FALSE
  )$RNA
  pb_ref <- AggregateExpression(
    ref_seurat,
    assays = "RNA",
    slot = "counts",
    features = common_genes,
    return.seurat = FALSE
  )$RNA
  keep_genes <- rowSums(pb_my) > 10 & rowSums(pb_ref) > 10
  pb_my <- pb_my[keep_genes, ]
  pb_ref <- pb_ref[keep_genes, ]
  normalize_pb <- function(pb_matrix) {
    cpm <- sweep(pb_matrix, 2, colSums(pb_matrix), FUN = "/") * 1e6
    log1p(cpm)
  }
  pb_my_norm <- normalize_pb(pb_my)
  pb_ref_norm <- normalize_pb(pb_ref)
  gene_vars_my <- apply(pb_my_norm, 1, var)
  gene_vars_ref <- apply(pb_ref_norm, 1, var)
  top_genes_my <- names(sort(gene_vars_my, decreasing = TRUE))[1:min(3000, length(gene_vars_my))]
  top_genes_ref <- names(sort(gene_vars_ref, decreasing = TRUE))[1:min(3000, length(gene_vars_ref))]
  features_to_use <- intersect(union(top_genes_my, top_genes_ref), common_genes)
  pb_my_norm_final <- pb_my_norm[features_to_use, ]
  pb_ref_norm_final <- pb_ref_norm[features_to_use, ]
  pb_my_norm_final <- t(scale(t(pb_my_norm_final)))
  pb_ref_norm_final <- t(scale(t(pb_ref_norm_final)))
  pb_my_norm_final[is.na(pb_my_norm_final)] <- 0
  pb_ref_norm_final[is.na(pb_ref_norm_final)] <- 0
  celltypes_my <- colnames(pb_my_norm_final)
  celltypes_ref <- colnames(pb_ref_norm_final)
  cor_matrix <- matrix(
    NA,
    nrow = length(celltypes_my),
    ncol = length(celltypes_ref),
    dimnames = list(celltypes_my, celltypes_ref)
  )
  pval_matrix <- cor_matrix
  for (i in celltypes_my) {
    for (j in celltypes_ref) {
      test <- cor.test(pb_my_norm_final[, i], pb_ref_norm_final[, j], method = "pearson")
      cor_matrix[i, j] <- test$estimate
      pval_matrix[i, j] <- test$p.value
    }
  }
  padj_matrix <- matrix(
    p.adjust(pval_matrix, method = "BH"),
    nrow = nrow(pval_matrix),
    dimnames = dimnames(pval_matrix)
  )
  cor_df <- as.data.frame(as.table(cor_matrix))
  colnames(cor_df) <- c("My_Data", "Ref_Data", "Correlation")
  padj_df <- as.data.frame(as.table(padj_matrix))
  colnames(padj_df) <- c("My_Data", "Ref_Data", "Padj")
  plot_data <- cor_df %>%
    left_join(padj_df, by = c("My_Data", "Ref_Data")) %>%
    mutate(Significance = ifelse(Padj < 0.05 & Correlation > 0.7, "*", ""))
  export_table <- plot_data %>%
    transmute(
      My_Atlas_CellType = My_Data,
      Reference_Atlas_CellType = Ref_Data,
      Pearson_Correlation = Correlation,
      Adjusted_P_Value = Padj,
      Significance_Mark = Significance
    ) %>%
    arrange(desc(Pearson_Correlation))
  source_csv <- paste0(output_base, "_source_data.csv")
  write.csv(export_table, source_csv, row.names = FALSE)
}
