suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})
input_rds <- "/Files/xiezhouwei/MS/Online_Sn/MS_Celltype_meta_final.rds"
old_selected_file <- "/data/work/RCTD_26.6/RCTD_RefAtlas_Downsample_noNEU_filtered/RCTD_reference_noNEU_2000perType_bySample_exclude_AS_Ma24_MS549H_MS549T_selected_cell_ids.txt"
old_summary_file <- "/data/work/RCTD_26.6/RCTD_RefAtlas_Downsample_noNEU_filtered/RCTD_reference_noNEU_2000perType_bySample_exclude_AS_Ma24_MS549H_MS549T_sampling_summary.csv"
output_dir <- "/data/work/RCTD_26.6/RCTD_GM_Ref_Downsample_withNEU_filtered"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
out_rds <- file.path(output_dir, "reference.rds")
out_meta <- file.path(output_dir, "metadata.csv")
out_summary <- file.path(output_dir, "sampling_summary.csv")
out_diag <- file.path(output_dir, "diagnostic_long.csv")
out_selected <- file.path(output_dir, "selected_cell_ids.txt")
out_check <- file.path(output_dir, "retained_noNEU_check.csv")
out_session <- file.path(output_dir, "sessionInfo.txt")
celltype_key <- "celltype_final"
sample_key <- "Sample_id"
keep_celltypes <- c("AS", "IMM", "MG", "NEU", "OL", "OPC", "Vas")
old_celltypes <- setdiff(keep_celltypes, "NEU")
new_celltype <- "NEU"
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
    members <- rownames(sub)[sample_ids == sid]
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
message("[Read old noNEU selected cells] ", old_selected_file)
old_selected_cells <- readLines(old_selected_file, warn = FALSE)
old_selected_cells <- old_selected_cells[nzchar(old_selected_cells)]
if (length(old_selected_cells) == 0) {
  stop("Old selected cell ID file is empty.")
}
if (anyDuplicated(old_selected_cells) > 0) {
  stop("Duplicated cell IDs found in old selected file.")
}
missing_old <- setdiff(old_selected_cells, colnames(obj))
if (length(missing_old) > 0) {
  stop("Some old selected cells are not present in input object. Example: ", missing_old[1])
}
old_meta <- meta[old_selected_cells, , drop = FALSE]
old_types_found <- unique(as.character(old_meta[[celltype_key]]))
unexpected_old_types <- setdiff(old_types_found, old_celltypes)
if (length(unexpected_old_types) > 0) {
  stop(
    "Old selected file contains unexpected cell types: ",
    paste(unexpected_old_types, collapse = ", ")
  )
}
old_as_bad <- (
  as.character(old_meta[[celltype_key]]) == "AS" &
    as.character(old_meta[[sample_key]]) %in% as_exclude_sample_ids
)
if (any(old_as_bad)) {
  stop("Old selected cells contain AS cells from excluded sample IDs.")
}
old_counts <- table(factor(as.character(old_meta[[celltype_key]]), levels = old_celltypes))
message("[Retained old noNEU cell counts]")
print(old_counts)
if (any(as.integer(old_counts) == 0)) {
  stop("Old noNEU selected file is missing at least one expected cell type.")
}
qc_keep <- (
  as.numeric(meta$percent.mt) <= max_mt_pct &
    as.numeric(meta$nFeature_RNA) >= min_genes &
    as.numeric(meta$nCount_RNA) >= min_counts
)
neu_eligible <- as.character(meta[[celltype_key]]) == new_celltype & qc_keep
meta_neu_eligible <- meta[neu_eligible, , drop = FALSE]
message("[NEU eligible count after QC]")
print(nrow(meta_neu_eligible))
neu_seed <- seed + which(keep_celltypes == new_celltype) * 1009
neu_res <- sample_one_celltype(
  meta = meta_neu_eligible,
  ct = new_celltype,
  target_n = target_per_celltype,
  seed_i = neu_seed
)
neu_cells <- neu_res$cells
if (length(intersect(old_selected_cells, neu_cells)) > 0) {
  stop("NEU selected cells overlap with old selected cells.")
}
selected_cells <- c(old_selected_cells, neu_cells)
if (anyDuplicated(selected_cells) > 0) {
  stop("Duplicated selected cell IDs detected after adding NEU.")
}
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
for (ct in old_celltypes) {
  old_ids_ct <- rownames(old_meta)[as.character(old_meta[[celltype_key]]) == ct]
  new_ids_ct <- colnames(obj_ref)[as.character(obj_ref@meta.data[[celltype_key]]) == ct]
  if (!setequal(old_ids_ct, new_ids_ct)) {
    stop("Retained old cell IDs changed for cell type: ", ct)
  }
}
if (as.integer(final_counts[new_celltype]) == 0) {
  stop("NEU has zero cells in final object.")
}
if (ncol(obj_ref) != length(selected_cells)) {
  stop(
    "Final Seurat object has ", ncol(obj_ref),
    " cells, but selected_cells has ", length(selected_cells)
  )
}
if (file.exists(old_summary_file)) {
  old_summary <- read.csv(old_summary_file, stringsAsFactors = FALSE, check.names = FALSE)
  old_summary <- old_summary[old_summary$celltype %in% old_celltypes, , drop = FALSE]
  old_summary$sampling_mode <- paste0(
    old_summary$sampling_mode,
    "; retained_from_previous_noNEU_run"
  )
} else {
  old_summary <- data.frame(
    celltype = old_celltypes,
    n_available_before_sample_id_filter = NA_integer_,
    n_available = as.integer(old_counts),
    n_sampled = as.integer(old_counts),
    sampling_mode = "retained_from_previous_noNEU_selected_cell_ids",
    n_sample_ids_available = as.integer(sapply(old_celltypes, function(ct) {
      length(unique(as.character(old_meta[as.character(old_meta[[celltype_key]]) == ct, sample_key])))
    })),
    n_sample_ids_sampled = as.integer(sapply(old_celltypes, function(ct) {
      length(unique(as.character(old_meta[as.character(old_meta[[celltype_key]]) == ct, sample_key])))
    }))
  )
}
summary_df <- rbind(old_summary, neu_res$summary)
rownames(summary_df) <- NULL
retained_check <- data.frame(
  celltype = old_celltypes,
  n_old = as.integer(old_counts),
  n_final = as.integer(final_counts[old_celltypes]),
  unchanged = as.integer(old_counts) == as.integer(final_counts[old_celltypes])
)
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
message("[Save] retained noNEU check: ", out_check)
write.csv(retained_check, out_check, row.names = FALSE)
writeLines(capture.output(sessionInfo()), out_session)
message("[Done]")
message("Output RDS: ", out_rds)
