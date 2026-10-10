suppressPackageStartupMessages({
  library(hdWGCNA)
  library(Seurat)
  library(WGCNA)
  library(dplyr)
  library(tidyverse)
})
set.seed(12345)
enableWGCNAThreads(nThreads = 16)
input_rds <- "/data/work/Online_data/2_AS_anlysis/hdWGCNA_AS_ALLregion/CombinedMS_AS_allregion.rds"
outdir <- "/data/work/Online_data/2_AS_anlysis/hdWGCNA_AS_ALLregion/AS_allregion_hdWGCNA_to_ConstructNetwork"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
wgcna_name <- "AS_allRegion_WGCNA"
tom_name <- "AS_allRegion"
fraction <- 0.05
checkpoint_metacells <- file.path(outdir, "checkpoint_01_after_metacells_before_softpower.rds")
checkpoint_softpower <- file.path(outdir, "checkpoint_02_after_softpower.rds")
final_network_rds <- file.path(outdir, "checkpoint_03_after_ConstructNetwork.rds")
seurat_as <- readRDS(input_rds)
if (!inherits(seurat_as, "Seurat")) {
  stop("Input object is not a Seurat object.")
}
seurat_as <- SeuratObject::UpdateSeuratObject(seurat_as)
DefaultAssay(seurat_as) <- "RNA"
required_cols <- c(
  "celltype_final", "Sample_id", "donor_id", "Study", "matter",
  "lesion_type", "diagnosis", "sex", "age_at_death",
  "nCount_RNA", "nFeature_RNA", "percent.mt"
)
missing_cols <- setdiff(required_cols, colnames(seurat_as@meta.data))
if (length(missing_cols) > 0) {
  stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "))
}
if (any(seurat_as$celltype_final != "AS")) {
  seurat_as <- subset(seurat_as, subset = celltype_final == "AS")
}
reduction_use <- if ("harmony" %in% names(seurat_as@reductions)) {
  "harmony"
} else if ("pca" %in% names(seurat_as@reductions)) {
  "pca"
} else {
  stop("No harmony or pca reduction found.")
}
seurat_as <- SetupForWGCNA(
  seurat_as,
  gene_select = "fraction",
  fraction = fraction,
  wgcna_name = wgcna_name
)
seurat_as <- MetacellsByGroups(
  seurat_obj = seurat_as,
  group.by = c("celltype_final", "Sample_id"),
  reduction = reduction_use,
  k = 25,
  max_shared = 10,
  ident.group = "celltype_final",
  min_cells = 50,
  wgcna_name = wgcna_name
)
seurat_as <- NormalizeMetacells(seurat_as, wgcna_name = wgcna_name)
saveRDS(seurat_as, checkpoint_metacells)
seurat_as <- SetDatExpr(
  seurat_as,
  group_name = "AS",
  group.by = "celltype_final",
  assay = "RNA",
  layer = "data",
  wgcna_name = wgcna_name
)
datExpr <- GetDatExpr(
  seurat_as,
  wgcna_name = wgcna_name
)
dat_dim <- data.frame(
  n_metacells = nrow(datExpr),
  n_genes_for_WGCNA = ncol(datExpr)
)
print(dat_dim)
write_csv(
  dat_dim,
  file.path(outdir, "02_datExpr_dimension_before_TestSoftPowers.csv")
)
cat("\nFirst few genes in datExpr:\n")
print(head(colnames(datExpr)))
cat("\nFirst few metacells in datExpr:\n")
print(head(rownames(datExpr)))
cat("\nRunning TestSoftPowers...\n")
cat("Start time:", as.character(Sys.time()), "\n")
seurat_as <- TestSoftPowers(
  seurat_as,
  networkType = "signed",
  wgcna_name = wgcna_name
)
cat("Finished TestSoftPowers.\n")
cat("End time:", as.character(Sys.time()), "\n")
power_table <- GetPowerTable(
  seurat_as,
  wgcna_name = wgcna_name
)
write_csv(
  power_table,
  file.path(outdir, "02_soft_power_table.csv")
)
print(power_table)
saveRDS(
  seurat_as,
  checkpoint_softpower
)
tom_outdir <- file.path(outdir, "TOM")
dir.create(tom_outdir, recursive = TRUE, showWarnings = FALSE)
wgcna_name <- "AS_allRegion_WGCNA"
tom_name <- "AS_allRegion"
construct_args <- names(formals(ConstructNetwork))
if ("power" %in% construct_args) {
  cat("ConstructNetwork supports `power`; using power = 7 explicitly.\n")
  seurat_as <- ConstructNetwork(
    seurat_as,
    tom_name = tom_name,
    tom_outdir = tom_outdir,
    networkType = "signed",
    power = 7,
    wgcna_name = wgcna_name,
    overwrite_tom = TRUE
  )
} else {
  cat("ConstructNetwork does not support `power`; using power stored in hdWGCNA object.\n")
  seurat_as <- ConstructNetwork(
    seurat_as,
    tom_name = tom_name,
    tom_outdir = tom_outdir,
    networkType = "signed",
    wgcna_name = wgcna_name,
    overwrite_tom = TRUE
  )
}
saveRDS(
  seurat_as,
  final_network_rds
)
