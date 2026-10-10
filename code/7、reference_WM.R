suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})
input_rds <- "/Files/xiezhouwei/MS/Online_Sn/MS_Celltype_meta_final.rds"
output_dir <- "/data/work/RCTD_26.6"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_prefix <- "RCTD_reference_noNEU_2000perType_bySample_exclude_AS_Ma24_MS549H_MS549T"
out_rds <- file.path(output_dir, paste0(output_prefix, ".rds"))
out_meta <- file.path(output_dir, paste0(output_prefix, "_metadata.csv"))
out_summary <- file.path(output_dir, paste0(output_prefix, "_sampling_summary.csv"))
out_diag <- file.path(output_dir, paste0(output_prefix, "_diagnostic_long.csv"))
out_selected <- file.path(output_dir, paste0(output_prefix, "_selected_cell_ids.txt"))
out_excluded_as <- file.path(output_dir, paste0(output_prefix, "_AS_excluded_cells_metadata.csv"))
out_session <- file.path(output_dir, paste0(output_prefix, "_sessionInfo.txt"))
celltype_key <- "celltype_final"
sample_key <- "Sample_id"
keep_celltypes <- c("AS", "IMM", "MG", "OL", "OPC", "Vas")
as_exclude_sample_ids <- c("Ma24_MS549H", "Ma24_MS549T")
target_per_celltype <- 2000
seed <- 20260626
max_mt_pct <- 10
min_genes <- 500
min_counts <- 1000
diagnostic_columns <- c(
  "Sample_id", "diagnosis", "lesion_type",
  "matter", "sex", "age_at_death"
)
get_counts <- function(obj, assay = DefaultAssay(obj)) {
  x <- tryCatch(
    GetAssayData(obj, assay = assay, layer = "counts"),
    error = function(e) NULL
  )
  if (is.null(x)) {
    x <- tryCatch(
      GetAssayData(obj, assay = assay, slot = "counts"),
      error = function(e) NULL
    )
  }
  if (is.null(x) || nrow(x) == 0 || ncol(x) == 0) {
    stop("No valid counts matrix found in assay: ", assay)
  }
  x
}
largest_remainder_allocation <- function(group_counts, target_n) {
  if (sum(group_counts) <= target_n) {
    return(group_counts)
  }
  expected <- group_counts / sum(group_counts) * target_n
  alloc <- floor(expected)
  remainder <- target_n - sum(alloc)
  if (remainder > 0) {
    frac_order <- names(sort(expected - alloc, decreasing = TRUE))
    for (sid in frac_order) {
      if (remainder == 0) break
      if (alloc[[sid]] < group_counts[[sid]]) {
        alloc[[sid]] <- alloc[[sid]] + 1L
        remainder <- remainder - 1L
      }
    }
  }
  alloc
}
sample_one_celltype <- function(meta, ct, target_n, seed_i) {
  sub <- meta[as.character(meta[[celltype_key]]) == ct, , drop = FALSE]
  n_available_before_sample_id_filter <- nrow(sub)
  if (nrow(sub) == 0) {
    stop("No cells available for cell type: ", ct)
  }
  sample_ids <- as.character(sub[[sample_key]])
  missing_sample_id <- is.na(sample_ids) | sample_ids == "" | sample_ids == "NA"
  if (any(missing_sample_id)) {
    message(
      "[Sample_id filter] ", ct, ": excluding ",
      sum(missing_sample_id),
      " cells with missing Sample_id before proportional sampling."
    )
    sub <- sub[!missing_sample_id, , drop = FALSE]
    sample_ids <- as.character(sub[[sample_key]])
  }
  n_available <- nrow(sub)
  if (n_available == 0) {
    stop("No cells with valid Sample_id available for cell type: ", ct)
  }
  if (n_available <= target_n) {
    warning(
      ct, " has only ", n_available,
      " QC-passed cells with valid Sample_id; keeping all available cells."
    )
    selected <- rownames(sub)
    return(list(
      cells = selected,
      summary = data.frame(
        celltype = ct,
        n_available_before_sample_id_filter = n_available_before_sample_id_filter,
        n_available = n_available,
        n_sampled = length(selected),
        sampling_mode = "all_available_QC_passed_valid_Sample_id",
        n_sample_ids_available = length(unique(sample_ids)),
        n_sample_ids_sampled = length(unique(as.character(sub[[sample_key]])))
      )
    ))
  }
  group_counts <- table(sample_ids)
  alloc <- largest_remainder_allocation(group_counts, target_n)
  set.seed(seed_i)
  selected <- unlist(lapply(names(alloc), function(sid) {
    n_take <- as.integer(alloc[[sid]])
    if (n_take <= 0) return(character(0))
    members <- rownames(sub)[!is.na(sample_ids) & sample_ids == sid]
    if (length(members) < n_take) {
      stop(
        "Internal allocation error for ", ct, " / ", sid,
        ": requested ", n_take, " but only ", length(members), " available."
      )
    }
    sample(members, n_take, replace = FALSE)
  }), use.names = FALSE)
  selected <- selected[!is.na(selected)]
  if (length(selected) != target_n) {
    stop(
      "Sampling failed for ", ct,
      ": selected ", length(selected),
      " cells, expected ", target_n
    )
  }
  sampled_sub <- sub[selected, , drop = FALSE]
  list(
    cells = selected,
    summary = data.frame(
      celltype = ct,
      n_available_before_sample_id_filter = n_available_before_sample_id_filter,
      n_available = n_available,
      n_sampled = length(selected),
      sampling_mode = "proportional_by_Sample_id_largest_remainder",
      n_sample_ids_available = length(unique(sample_ids)),
      n_sample_ids_sampled = length(unique(as.character(sampled_sub[[sample_key]])))
    )
  )
}
make_diagnostic_long <- function(meta) {
  present_cols <- intersect(diagnostic_columns, colnames(meta))
  if (length(present_cols) == 0) {
    return(data.frame(
      celltype_final = character(),
      diagnostic_column = character(),
      level = character(),
      n = integer()
    ))
  }
  rows <- lapply(present_cols, function(col) {
    tab <- as.data.frame(table(
      celltype_final = as.character(meta[[celltype_key]]),
      level = as.character(meta[[col]]),
      useNA = "ifany"
    ))
    tab <- tab[tab$Freq > 0, , drop = FALSE]
    data.frame(
      celltype_final = tab$celltype_final,
      diagnostic_column = col,
      level = as.character(tab$level),
      n = as.integer(tab$Freq)
    )
  })
  do.call(rbind, rows)
}
set.seed(seed)
message("[Load] Reading RDS: ", input_rds)
obj <- readRDS(input_rds)
if (!inherits(obj, "Seurat")) {
  stop("Input is not a Seurat object. Class: ", paste(class(obj), collapse = ", "))
}
meta <- obj@meta.data
required_cols <- c(celltype_key, sample_key)
missing_cols <- setdiff(required_cols, colnames(meta))
if (length(missing_cols) > 0) {
  stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "))
}
assay <- DefaultAssay(obj)
counts <- get_counts(obj, assay)
message("[Input] assay: ", assay)
message("[Input] cells: ", ncol(obj), "; genes: ", nrow(obj))
if (!identical(rownames(meta), colnames(obj))) {
  stop("Metadata rownames are not identical to Seurat object colnames.")
}
if (!"nFeature_RNA" %in% colnames(meta)) {
  meta$nFeature_RNA <- Matrix::colSums(counts > 0)
  obj$nFeature_RNA <- meta$nFeature_RNA
}
if (!"nCount_RNA" %in% colnames(meta)) {
  meta$nCount_RNA <- Matrix::colSums(counts)
  obj$nCount_RNA <- meta$nCount_RNA
}
if (!"percent.mt" %in% colnames(meta)) {
  mt_genes <- grep("^MT-|^mt-", rownames(counts), value = TRUE)
  if (length(mt_genes) > 0) {
    total_counts <- Matrix::colSums(counts)
    mt_counts <- Matrix::colSums(counts[mt_genes, , drop = FALSE])
    meta$percent.mt <- ifelse(total_counts > 0, mt_counts / total_counts * 100, 0)
  } else {
    meta$percent.mt <- 0
  }
  obj$percent.mt <- meta$percent.mt
}
message("[Input celltype counts]")
print(table(meta[[celltype_key]], useNA = "ifany"))
base_keep <- as.character(meta[[celltype_key]]) %in% keep_celltypes
qc_keep <- (
  as.numeric(meta$percent.mt) <= max_mt_pct &
    as.numeric(meta$nFeature_RNA) >= min_genes &
    as.numeric(meta$nCount_RNA) >= min_counts
)
as_exclude <- (
  as.character(meta[[celltype_key]]) == "AS" &
    as.character(meta[[sample_key]]) %in% as_exclude_sample_ids
)
excluded_as_meta <- meta[as_exclude, , drop = FALSE]
if (nrow(excluded_as_meta) > 0) {
  write.csv(
    data.frame(cell_id = rownames(excluded_as_meta), excluded_as_meta, check.names = FALSE),
    out_excluded_as,
    row.names = FALSE
  )
}
eligible <- base_keep & qc_keep & !as_exclude
meta_eligible <- meta[eligible, , drop = FALSE]
message("[Eligible counts after NEU removal, QC, and AS exclusion]")
print(table(factor(as.character(meta_eligible[[celltype_key]]), levels = keep_celltypes)))
message("[Missing Sample_id counts among eligible cells]")
print(tapply(
  is.na(meta_eligible[[sample_key]]) |
    as.character(meta_eligible[[sample_key]]) == "" |
    as.character(meta_eligible[[sample_key]]) == "NA",
  factor(as.character(meta_eligible[[celltype_key]]), levels = keep_celltypes),
  sum
))
selected_cells <- character(0)
summary_list <- list()
for (i in seq_along(keep_celltypes)) {
  ct <- keep_celltypes[i]
  res <- sample_one_celltype(
    meta = meta_eligible,
    ct = ct,
    target_n = target_per_celltype,
    seed_i = seed + i * 1009
  )
  selected_cells <- c(selected_cells, res$cells)
  summary_list[[ct]] <- res$summary
}
selected_cells <- selected_cells[!is.na(selected_cells)]
if (length(selected_cells) == 0) {
  stop("No selected cells.")
}
if (anyDuplicated(selected_cells) > 0) {
  stop("Duplicated selected cell IDs detected.")
}
if (any(is.na(selected_cells))) {
  stop("NA selected cell IDs detected.")
}
summary_df <- do.call(rbind, summary_list)
rownames(summary_df) <- NULL
message("[Sampling summary]")
print(summary_df)
message("[Subset] Creating final Seurat reference object")
obj_ref <- subset(obj, cells = selected_cells)
obj_ref[[celltype_key]] <- factor(
  as.character(obj_ref@meta.data[[celltype_key]]),
  levels = keep_celltypes
)
final_counts <- table(
  factor(as.character(obj_ref@meta.data[[celltype_key]]), levels = keep_celltypes)
)
message("[Final selected counts]")
print(final_counts)
if (any(as.integer(final_counts) == 0)) {
  stop("At least one target cell type has zero cells in final object.")
}
if (ncol(obj_ref) != length(selected_cells)) {
  stop(
    "Final Seurat object has ", ncol(obj_ref),
    " cells, but selected_cells has ", length(selected_cells)
  )
}
if (any(is.na(colnames(obj_ref)))) {
  stop("Final Seurat object contains NA cell names.")
}
as_bad <- (
  as.character(obj_ref@meta.data[[celltype_key]]) == "AS" &
    as.character(obj_ref@meta.data[[sample_key]]) %in% as_exclude_sample_ids
)
if (any(as_bad)) {
  stop("AS excluded sample IDs are still present in final object.")
}
meta_ref <- data.frame(
  cell_id = colnames(obj_ref),
  obj_ref@meta.data,
  check.names = FALSE
)
diag_long <- make_diagnostic_long(obj_ref@meta.data)
message("[Save] RDS: ", out_rds)
saveRDS(obj_ref, out_rds, compress = "xz")
message("[Save] metadata: ", out_meta)
write.csv(meta_ref, out_meta, row.names = FALSE)
message("[Save] sampling summary: ", out_summary)
write.csv(summary_df, out_summary, row.names = FALSE)
message("[Save] diagnostic long: ", out_diag)
write.csv(diag_long, out_diag, row.names = FALSE)
message("[Save] selected cell ids: ", out_selected)
writeLines(selected_cells, out_selected)
writeLines(capture.output(sessionInfo()), out_session)
message("[Done]")
message("Output RDS: ", out_rds)
