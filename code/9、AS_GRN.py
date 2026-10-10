import os
import gc
import pickle
import warnings
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse import csr_matrix
import stereo as st
from stereo.io.reader import read_h5ms
from stereo.io.writer import write_h5ms
from stereo.core.ms_data import MSData
from stereo.core.ms_pipeline import slice_generator as sg
from ctxcore.rnkdb import (
    FeatherRankingDatabase as RankingDatabase
)
ms_data = st.io.read_h5ms('/data/work/Stereopy_CCC_NEW/WM_noNEU/MS_microenv_with_CCC.h5ms')
import os, gc, warnings
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse import csr_matrix
from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
GRN_DIR = "/data/work/GRN/3D_GRN_human"
TF_FILE = os.path.join(GRN_DIR, "allTFs_hg38_in_data.txt")
DATABASE = os.path.join(GRN_DIR, "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")
MOTIF_ANNO = os.path.join(GRN_DIR, "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl")
OUT_DIR = "/data/work/GRN/AS_multi_GRN"
PREPROCESSED_H5MS = os.path.join(OUT_DIR, "AS_multi_rawcounts_selected_genes.h5ms")
FINAL_H5MS = os.path.join(OUT_DIR, "AS_multi_GRN_result.h5ms")
os.makedirs(OUT_DIR, exist_ok=True)
for path in [TF_FILE, DATABASE, MOTIF_ANNO]:
    if not os.path.exists(path):
        raise FileNotFoundError(path)
TARGET_CELLTYPE = "AS"
CELLTYPE_COL = "celltype_final"
MAX_AS_PER_SAMPLE = None
N_BASE_GENES = 5000
MIN_GENE_CELLS_TOTAL = 20
MIN_GENE_SAMPLES = 2
MIN_TF_CELLS_TOTAL = 10
MIN_TF_SAMPLES = 2
N_WORKERS = 8
RANDOM_SEED = 20260716
def as_str_array(values):
    return np.asarray(values).astype(str)
def looks_numeric_placeholder(names, n_check=200):
    names = as_str_array(names)
    if len(names) == 0:
        return False
    names = names[:min(n_check, len(names))]
    return np.mean([x.isdigit() for x in names]) >= 0.95
def get_effective_raw_axes(data, sample_name):
    if data.raw is None:
        raise RuntimeError(f"{sample_name}: raw不存在。")
    current_cells, current_genes = as_str_array(data.cell_names), as_str_array(data.gene_names)
    raw_cells, raw_genes = as_str_array(data.raw.cell_names), as_str_array(data.raw.gene_names)
    if data.raw.exp_matrix.shape != (len(current_cells), len(current_genes)):
        raise RuntimeError(f"{sample_name}: raw矩阵维度与当前对象不一致。")
    if np.array_equal(raw_genes, current_genes):
        effective_genes, gene_status = raw_genes, "matched"
    elif looks_numeric_placeholder(raw_genes) and len(raw_genes) == len(current_genes):
        effective_genes, gene_status = current_genes, "recovered_by_current_order"
    elif set(raw_genes) == set(current_genes):
        effective_genes, gene_status = raw_genes, "raw_order_retained"
    else:
        raise RuntimeError(f"{sample_name}: 无法确定raw基因对应关系。")
    if np.array_equal(raw_cells, current_cells):
        effective_cells, cell_status = raw_cells, "matched"
    elif looks_numeric_placeholder(raw_cells) and len(raw_cells) == len(current_cells):
        effective_cells, cell_status = current_cells, "recovered_by_current_order"
    elif set(raw_cells) == set(current_cells):
        effective_cells, cell_status = raw_cells, "raw_order_retained"
    else:
        raise RuntimeError(f"{sample_name}: 无法确定raw细胞对应关系。")
    return {"raw_cells": effective_cells, "raw_genes": effective_genes, "gene_status": gene_status, "cell_status": cell_status}
def sparse_subset(matrix, row_indices, col_indices):
    row_indices, col_indices = np.asarray(row_indices, dtype=int), np.asarray(col_indices, dtype=int)
    try:
        subset = matrix[row_indices, :][:, col_indices]
    except Exception:
        subset = np.asarray(matrix)[np.ix_(row_indices, col_indices)]
    return subset.tocsr() if sparse.issparse(subset) else csr_matrix(np.asarray(subset))
with open(TF_FILE, "r", encoding="utf-8") as f:
    all_tfs = list(dict.fromkeys(line.strip() for line in f if line.strip()))
ranking_db = RankingDatabase(fname=DATABASE, name="hg38_10kb")
db_genes = set(map(str, ranking_db.genes))
print("Samples:", list(ms_data.names))
print("Sample number:", ms_data.num_slice)
print("TF file:", len(all_tfs))
print("Ranking database genes:", len(db_genes))
preflight_rows, sample_gene_sets = [], []
for sample_idx, sample_name in enumerate(ms_data.names):
    data = ms_data[sample_name]
    meta = data.cells.to_df().copy()
    meta.index = as_str_array(data.cell_names)
    if CELLTYPE_COL not in meta.columns:
        raise KeyError(f"{sample_name}: 不存在元数据列 {CELLTYPE_COL}")
    genes = as_str_array(data.gene_names)
    if len(np.unique(genes)) != len(genes):
        raise RuntimeError(f"{sample_name}: gene_names存在重复。")
    raw_axes = get_effective_raw_axes(data, sample_name)
    gene_set = set(genes)
    n_as = int((meta[CELLTYPE_COL].astype(str) == TARGET_CELLTYPE).sum())
    tf_matches = sum(tf in gene_set for tf in all_tfs)
    db_matches = sum(gene in db_genes for gene in genes)
    preflight_rows.append({
        "sample_index": sample_idx, "sample": sample_name, "n_cells": data.n_cells, "n_genes": data.n_genes,
        "n_AS": n_as, "raw_gene_status": raw_axes["gene_status"], "raw_cell_status": raw_axes["cell_status"],
        "tf_matches": tf_matches, "db_gene_matches": db_matches, "db_match_fraction": db_matches / len(genes)
    })
    sample_gene_sets.append(gene_set)
preflight = pd.DataFrame(preflight_rows)
print("\nPreflight summary:")
print(preflight.to_string(index=False))
preflight.to_csv(os.path.join(OUT_DIR, "01_AS_multi_sample_preflight.csv"), index=False)
if (preflight["n_AS"] == 0).any():
    raise RuntimeError("以下样本没有AS：" + ", ".join(preflight.loc[preflight["n_AS"] == 0, "sample"]))
common_gene_set = set.intersection(*sample_gene_sets)
first_genes = as_str_array(ms_data[ms_data.names[0]].gene_names)
common_genes = [gene for gene in first_genes if gene in common_gene_set]
common_tf_names = [tf for tf in all_tfs if tf in common_gene_set]
common_db_genes = [gene for gene in common_genes if gene in db_genes]
print("\nCross-sample matching:")
print("Common genes:", len(common_genes))
print("Common TFs:", len(common_tf_names))
print("Common genes in ranking DB:", len(common_db_genes), f"({len(common_db_genes) / len(common_genes):.2%})")
print("First common TFs:", common_tf_names[:20])
if len(common_tf_names) < 100:
    raise RuntimeError("跨样本共有TF少于100。")
if len(common_db_genes) / len(common_genes) < 0.7:
    raise RuntimeError("跨样本共有基因与ranking数据库匹配率低于70%。")
with open(os.path.join(OUT_DIR, "02_common_genes_all_samples.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(common_genes) + "\n")
gc.collect()
import gc, os, warnings
import numpy as np
import pandas as pd
from scipy import sparse
from stereo.core.ms_data import MSData
from stereo.core.ms_pipeline import slice_generator as sg
rng = np.random.default_rng(RANDOM_SEED)
n_common_genes = len(common_genes)
total_sum = np.zeros(n_common_genes, dtype=np.float64)
total_sum_sq = np.zeros(n_common_genes, dtype=np.float64)
total_detected_cells = np.zeros(n_common_genes, dtype=np.int64)
detected_sample_number = np.zeros(n_common_genes, dtype=np.int16)
sample_selection, sample_count_rows = {}, []
total_as_cells = 0
for sample_idx, sample_name in enumerate(ms_data.names):
    data = ms_data[sample_name]
    meta = data.cells.to_df().copy()
    meta.index = as_str_array(data.cell_names)
    as_cells = meta.index[meta[CELLTYPE_COL].astype(str) == TARGET_CELLTYPE].tolist()
    original_n_as = len(as_cells)
    if MAX_AS_PER_SAMPLE is not None and len(as_cells) > MAX_AS_PER_SAMPLE:
        as_cells = rng.choice(as_cells, size=MAX_AS_PER_SAMPLE, replace=False).tolist()
    if len(as_cells) == 0:
        raise RuntimeError(f"{sample_name}: 没有可用AS。")
    sample_selection[sample_name] = as_cells
    raw_axes = get_effective_raw_axes(data, sample_name)
    raw_cell_index = pd.Index(raw_axes["raw_cells"])
    raw_gene_index = pd.Index(raw_axes["raw_genes"])
    row_idx = raw_cell_index.get_indexer(as_cells)
    col_idx = raw_gene_index.get_indexer(common_genes)
    if np.any(row_idx < 0):
        missing = np.asarray(as_cells)[row_idx < 0][:10]
        raise RuntimeError(f"{sample_name}: AS细胞未在raw中找到：{missing.tolist()}")
    if np.any(col_idx < 0):
        missing = np.asarray(common_genes)[col_idx < 0][:10]
        raise RuntimeError(f"{sample_name}: 共有基因未在raw中找到：{missing.tolist()}")
    raw_as = sparse_subset(data.raw.exp_matrix, row_idx, col_idx)
    expected_shape = (len(as_cells), len(common_genes))
    if raw_as.shape != expected_shape:
        raise RuntimeError(f"{sample_name}: raw AS矩阵维度{raw_as.shape}，预期{expected_shape}。")
    detected = np.asarray((raw_as > 0).sum(axis=0)).ravel()
    gene_sum = np.asarray(raw_as.sum(axis=0)).ravel()
    gene_sum_sq = np.asarray(raw_as.power(2).sum(axis=0)).ravel()
    total_detected_cells += detected.astype(np.int64)
    detected_sample_number += (detected > 0).astype(np.int16)
    total_sum += gene_sum
    total_sum_sq += gene_sum_sq
    total_as_cells += raw_as.shape[0]
    sample_count_rows.append({
        "sample_index": sample_idx, "sample": sample_name, "original_AS_cells": original_n_as,
        "selected_AS_cells": len(as_cells), "common_genes": raw_as.shape[1],
        "raw_nonzero": int(raw_as.nnz), "raw_total_counts": float(raw_as.sum())
    })
    print(sample_name, "AS:", len(as_cells), "raw matrix:", raw_as.shape, "nonzero:", raw_as.nnz)
    del raw_as
    gc.collect()
sample_counts = pd.DataFrame(sample_count_rows)
sample_counts.to_csv(os.path.join(OUT_DIR, "03_AS_cell_numbers_by_sample.csv"), index=False)
print("\nAS cell summary:")
print(sample_counts.to_string(index=False))
print("\nTotal AS cells:", total_as_cells)
if total_as_cells != sample_counts["selected_AS_cells"].sum():
    raise RuntimeError("累计AS细胞数与各样本AS细胞数之和不一致。")
imbalance_ratio = sample_counts["selected_AS_cells"].max() / sample_counts["selected_AS_cells"].min()
print("Maximum/minimum AS ratio:", round(imbalance_ratio, 3))
if imbalance_ratio > 3:
    warnings.warn("样本间AS细胞数相差超过3倍，较大样本可能对联合GRN贡献更多。")
global_mean = total_sum / total_as_cells
global_variance = (total_sum_sq - total_sum ** 2 / total_as_cells) / (total_as_cells - 1)
global_variance = np.maximum(global_variance, 0)
global_dispersion = global_variance / (global_mean + 1e-8)
tf_set = set(all_tfs)
gene_stats = pd.DataFrame({
    "detected_cells_total": total_detected_cells,
    "detected_samples": detected_sample_number,
    "sum_counts": total_sum,
    "mean_counts": global_mean,
    "variance": global_variance,
    "dispersion": global_dispersion,
    "in_ranking_database": [g in db_genes for g in common_genes],
    "is_TF": [g in tf_set for g in common_genes]
}, index=common_genes)
gene_stats.index.name = "gene"
eligible_genes = gene_stats.index[
    (gene_stats["detected_cells_total"] >= MIN_GENE_CELLS_TOTAL) &
    (gene_stats["detected_samples"] >= MIN_GENE_SAMPLES) &
    (gene_stats["variance"] > 0) &
    gene_stats["in_ranking_database"]
].tolist()
selected_tfs = [
    tf for tf in common_tf_names
    if gene_stats.loc[tf, "detected_cells_total"] >= MIN_TF_CELLS_TOTAL
    and gene_stats.loc[tf, "detected_samples"] >= MIN_TF_SAMPLES
]
selected_genes = list(dict.fromkeys(eligible_genes + selected_tfs))
selected_gene_set = set(selected_genes)
db_matches = [g for g in selected_genes if g in db_genes]
db_fraction = len(db_matches) / len(selected_genes)
print("\nGene selection summary:")
print("Common genes:", len(common_genes))
print("Eligible target genes:", len(eligible_genes))
print("Stable TFs:", len(selected_tfs))
print("Final selected genes:", len(selected_genes))
print("Genes in ranking DB:", len(db_matches), f"({db_fraction:.2%})")
if len(eligible_genes) < 1000:
    raise RuntimeError("满足条件的候选靶基因少于1000。")
if len(selected_tfs) < 100:
    raise RuntimeError("稳定表达TF少于100。")
if len(selected_genes) != len(selected_gene_set):
    raise RuntimeError("最终基因列表存在重复。")
if db_fraction < 0.70:
    raise RuntimeError("最终基因集与ranking database匹配率低于70%。")
gene_stats.to_csv(os.path.join(OUT_DIR, "04_AS_multi_gene_statistics.csv"))
with open(os.path.join(OUT_DIR, "05_AS_multi_selected_genes.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(selected_genes) + "\n")
with open(os.path.join(OUT_DIR, "06_AS_multi_selected_tfs.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(selected_tfs) + "\n")
as_data_list, as_sample_names, build_rows = [], [], []
for sample_idx, sample_name in enumerate(ms_data.names):
    data = ms_data[sample_name]
    selected_cells = sample_selection[sample_name]
    as_data = data.sub_by_name(
        cell_name=selected_cells,
        gene_name=selected_genes,
        filter_raw=False,
        copy=True
    )
    current_cells = as_str_array(as_data.cell_names)
    current_genes = as_str_array(as_data.gene_names)
    if as_data.n_cells != len(selected_cells):
        raise RuntimeError(f"{sample_name}: 子对象n_cells仍是旧值。")
    if as_data.n_genes != len(selected_genes):
        raise RuntimeError(f"{sample_name}: 子对象n_genes仍是旧值。")
    if as_data.exp_matrix.shape != (as_data.n_cells, as_data.n_genes):
        raise RuntimeError(f"{sample_name}: subset后exp_matrix维度与对象维度不一致。")
    if set(current_cells) != set(selected_cells):
        raise RuntimeError(f"{sample_name}: subset后细胞ID集合不一致。")
    if set(current_genes) != selected_gene_set:
        missing = sorted(selected_gene_set - set(current_genes))[:10]
        extra = sorted(set(current_genes) - selected_gene_set)[:10]
        raise RuntimeError(f"{sample_name}: subset后基因集合不一致；missing={missing}, extra={extra}")
    raw_axes = get_effective_raw_axes(data, sample_name)
    raw_cell_index = pd.Index(raw_axes["raw_cells"])
    raw_gene_index = pd.Index(raw_axes["raw_genes"])
    row_idx = raw_cell_index.get_indexer(current_cells)
    col_idx = raw_gene_index.get_indexer(current_genes)
    if np.any(row_idx < 0) or np.any(col_idx < 0):
        raise RuntimeError(f"{sample_name}: 子对象细胞或基因无法重新匹配raw。")
    raw_counts = sparse_subset(data.raw.exp_matrix, row_idx, col_idx)
    expected_shape = (len(current_cells), len(current_genes))
    if raw_counts.shape != expected_shape:
        raise RuntimeError(f"{sample_name}: raw counts维度{raw_counts.shape}，预期{expected_shape}。")
    as_data.exp_matrix = raw_counts.tocsr()
    if as_data.exp_matrix.shape != (as_data.n_cells, as_data.n_genes):
        raise RuntimeError(f"{sample_name}: 写入raw counts后矩阵维度与n_cells/n_genes不一致。")
    if len(as_data.cell_names) != as_data.exp_matrix.shape[0]:
        raise RuntimeError(f"{sample_name}: cell_names长度与矩阵行数不一致。")
    if len(as_data.gene_names) != as_data.exp_matrix.shape[1]:
        raise RuntimeError(f"{sample_name}: gene_names长度与矩阵列数不一致。")
    values = as_data.exp_matrix.data
    integer_fraction = np.mean(np.isclose(values, np.round(values), atol=1e-8)) if len(values) > 0 else 1.0
    if integer_fraction < 0.95:
        raise RuntimeError(f"{sample_name}: 提取后的表达矩阵不是原始整数counts。")
    as_data.cells["sample_name"] = np.repeat(sample_name, as_data.n_cells)
    as_data.cells["original_cell_id"] = current_cells
    as_data.cells[CELLTYPE_COL] = np.repeat(TARGET_CELLTYPE, as_data.n_cells)
    as_data_list.append(as_data)
    as_sample_names.append(sample_name)
    build_rows.append({
        "sample": sample_name, "n_cells": as_data.n_cells, "n_genes": as_data.n_genes,
        "matrix_rows": as_data.exp_matrix.shape[0], "matrix_columns": as_data.exp_matrix.shape[1],
        "nonzero": int(as_data.exp_matrix.nnz), "total_counts": float(as_data.exp_matrix.sum()),
        "integer_fraction": integer_fraction
    })
    print(sample_name, "final AS object:", as_data.n_cells, "x", as_data.n_genes)
    del raw_counts
    gc.collect()
build_summary = pd.DataFrame(build_rows)
build_summary.to_csv(os.path.join(OUT_DIR, "07_AS_object_build_summary.csv"), index=False)
print("\nFinal per-sample object check:")
print(build_summary.to_string(index=False))
if build_summary["n_genes"].nunique() != 1:
    raise RuntimeError("各样本AS对象的基因数不一致。")
if build_summary["matrix_columns"].nunique() != 1:
    raise RuntimeError("各样本表达矩阵列数不一致。")
if not np.all(build_summary["n_cells"] == build_summary["matrix_rows"]):
    raise RuntimeError("至少一个样本的n_cells与矩阵行数不一致。")
if not np.all(build_summary["n_genes"] == build_summary["matrix_columns"]):
    raise RuntimeError("至少一个样本的n_genes与矩阵列数不一致。")
as_ms_data = MSData(
    _data_list=as_data_list,
    _names=as_sample_names,
    _var_type="intersect",
    _relationship="other"
)
merged_as = as_ms_data.integrate()
expected_merged_cells = int(build_summary["n_cells"].sum())
expected_merged_genes = len(selected_genes)
expected_merged_shape = (expected_merged_cells, expected_merged_genes)
print("\nMerged AS object:")
print("n_cells:", merged_as.n_cells)
print("n_genes:", merged_as.n_genes)
print("exp_matrix shape:", merged_as.exp_matrix.shape)
print("expected shape:", expected_merged_shape)
if merged_as.n_cells != expected_merged_cells:
    raise RuntimeError("整合后的n_cells与各样本AS细胞数之和不一致。")
if merged_as.n_genes != expected_merged_genes:
    raise RuntimeError("整合后的n_genes与筛选基因数不一致。")
if merged_as.exp_matrix.shape != expected_merged_shape:
    raise RuntimeError("整合后的exp_matrix仍保留旧维度或整合不完整。")
if len(merged_as.cell_names) != expected_merged_cells:
    raise RuntimeError("整合后的cell_names长度与矩阵行数不一致。")
if len(merged_as.gene_names) != expected_merged_genes:
    raise RuntimeError("整合后的gene_names长度与矩阵列数不一致。")
if set(as_str_array(merged_as.gene_names)) != selected_gene_set:
    raise RuntimeError("整合后的基因集合与selected_genes不一致。")
print("\nAll matrix dimensions and axes are consistent.")
print("Final AS matrix:", merged_as.exp_matrix.shape)
print("Final TF number:", len(selected_tfs))
gc.collect()
import os, gc
import numpy as np
import pandas as pd
from scipy import sparse
from stereo.core.ms_pipeline import slice_generator as sg
RES_KEY = "AS_multi_regulatory_network_inference"
PREFIX = os.path.join(OUT_DIR, "AS_multi_GRN")
as_ms_data.to_integrate(res_key="sample_name", scope=sg[:], _from=sg[:], type="obs", item=["sample_name"] * as_ms_data.num_slice, cluster=False)
as_ms_data.to_integrate(res_key="original_cell_id", scope=sg[:], _from=sg[:], type="obs", item=["original_cell_id"] * as_ms_data.num_slice, cluster=False)
as_ms_data.to_integrate(res_key=CELLTYPE_COL, scope=sg[:], _from=sg[:], type="obs", item=[CELLTYPE_COL] * as_ms_data.num_slice, cluster=True)
merged_as = as_ms_data.merged_data
merged_meta = merged_as.cells.to_df().copy()
merged_meta.index = as_str_array(merged_as.cell_names)
merged_meta.index.name = "integrated_cell_id"
required_meta = ["batch", "sample_name", "original_cell_id", CELLTYPE_COL]
missing_meta = [x for x in required_meta if x not in merged_meta.columns]
if missing_meta:
    raise RuntimeError(f"整合对象缺少metadata：{missing_meta}")
if merged_meta["sample_name"].isna().any() or merged_meta["original_cell_id"].isna().any():
    raise RuntimeError("整合后部分细胞缺少sample_name或original_cell_id。")
if merged_meta.shape[0] != merged_as.n_cells:
    raise RuntimeError("metadata行数与整合矩阵细胞数不一致。")
if merged_meta.index.duplicated().any():
    raise RuntimeError("integrated_cell_id存在重复。")
cell_key = merged_meta["sample_name"].astype(str) + "||" + merged_meta["original_cell_id"].astype(str)
if cell_key.duplicated().any():
    raise RuntimeError("sample_name + original_cell_id组合存在重复。")
merged_meta[required_meta].to_csv(os.path.join(OUT_DIR, "08_AS_integrated_cell_id_map.csv"))
merged_genes = as_str_array(merged_as.gene_names)
merged_gene_set = set(merged_genes)
run_tfs = [tf for tf in selected_tfs if tf in merged_gene_set]
missing_tfs = sorted(set(selected_tfs) - merged_gene_set)
print("Merged matrix:", merged_as.exp_matrix.shape)
print("Metadata rows:", merged_meta.shape[0])
print("Selected TFs:", len(selected_tfs))
print("TFs passed to GRN:", len(run_tfs))
print("Missing TFs:", len(missing_tfs))
if merged_as.exp_matrix.shape != (merged_as.n_cells, merged_as.n_genes):
    raise RuntimeError("整合矩阵维度与n_cells/n_genes不一致。")
if len(run_tfs) != len(selected_tfs):
    raise RuntimeError(f"整合后有TF丢失：{missing_tfs[:20]}")
values = merged_as.exp_matrix.data if sparse.issparse(merged_as.exp_matrix) else np.asarray(merged_as.exp_matrix).ravel()
integer_fraction = np.mean(np.isclose(values, np.round(values), atol=1e-8)) if len(values) else 1.0
print("Integer fraction:", integer_fraction)
if integer_fraction < 0.95:
    raise RuntimeError("整合矩阵不再表现为原始整数counts。")
dense_gb = merged_as.n_cells * merged_as.n_genes * 8 / 1024**3
print(f"Estimated one float64 dense matrix: {dense_gb:.2f} GB")
print("\nStarting multi-sample AS-GRN...")
print("Samples:", as_ms_data.num_slice)
print("Cells:", merged_as.n_cells)
print("Genes:", merged_as.n_genes)
print("TFs:", len(run_tfs))
print("Workers:", N_WORKERS)
print("Result key:", RES_KEY)
import numpy as np
import pandas as pd
meta = merged_as.cells.to_df().copy()
batch = meta["batch"].astype(int).to_numpy()
old_position = np.asarray(merged_as.position, dtype=float)
new_position = old_position.copy()
batch_ids = sorted(np.unique(batch))
max_width = max(np.ptp(old_position[batch == b, 0]) for b in batch_ids)
max_height = max(np.ptp(old_position[batch == b, 1]) for b in batch_ids)
gap = max(max_width, max_height) * 10
cursor = 0.0
rows = []
for b in batch_ids:
    idx = batch == b
    x = old_position[idx, 0]
    y = old_position[idx, 1]
    shifted_x = x - x.min() + cursor
    shifted_y = y - y.min()
    new_position[idx, 0] = shifted_x
    new_position[idx, 1] = shifted_y
    rows.append({
        "batch": b,
        "n_cells": int(idx.sum()),
        "new_x_min": float(shifted_x.min()),
        "new_x_max": float(shifted_x.max()),
        "new_y_min": float(shifted_y.min()),
        "new_y_max": float(shifted_y.max())
    })
    cursor = shifted_x.max() + gap
merged_as.position = new_position
position_summary = pd.DataFrame(rows)
print(position_summary.to_string(index=False))
print("Gap:", gap)
for i in range(len(position_summary) - 1):
    current_max = position_summary.loc[i, "new_x_max"]
    next_min = position_summary.loc[i + 1, "new_x_min"]
    if next_min - current_max < gap * 0.9:
        raise RuntimeError("样本间距不足。")
print("All samples are spatially separated.")
HOTSPOT_KEY = "AS_multi_hotspot_regulatory_network_inference"
HOTSPOT_PREFIX = os.path.join(OUT_DIR, "AS_multi_hotspot_GRN")
print("\nStarting multi-sample AS Hotspot-GRN...")
print("Samples:", as_ms_data.num_slice)
print("Cells:", merged_as.n_cells)
print("Genes:", merged_as.n_genes)
print("TFs:", len(run_tfs))
print("Workers:", N_WORKERS)
print("Result key:", HOTSPOT_KEY)
as_ms_data.tl.regulatory_network_inference(
    scope=sg[:], mode="integrate",
    database=DATABASE, motif_anno=MOTIF_ANNO,
    tfs=run_tfs, target_genes=None,
    auc_threshold=0.05, num_workers=N_WORKERS,
    res_key=HOTSPOT_KEY, seed=RANDOM_SEED,
    cache=False, save_regulons=True,
    save_loom=False, fn_prefix=HOTSPOT_PREFIX,
    method="hotspot", ThreeD_slice=False,
    hotspot_kwargs={
        "model": "bernoulli",
        "n_neighbors": 30,
        "fdr_threshold": 0.05,
        "weighted_graph": False
    },
    use_raw=False
)
if HOTSPOT_KEY not in as_ms_data.merged_data.tl.result:
    raise RuntimeError(f"运行结束后未找到结果键：{HOTSPOT_KEY}")
result = as_ms_data.merged_data.tl.result[HOTSPOT_KEY]
required_keys = ["adjacencies", "auc_matrix", "regulons"]
missing_keys = [key for key in required_keys if key not in result]
if missing_keys:
    raise RuntimeError(f"GRN结果缺少：{missing_keys}")
adj = result["adjacencies"]
auc = result["auc_matrix"]
regulons = result["regulons"]
print("\nHotspot-GRN completed.")
print("Result keys:", list(result.keys()))
print("Adjacencies:", adj.shape)
print("AUC matrix:", auc.shape)
print("Regulons:", len(regulons))
if auc.shape[0] != merged_as.n_cells:
    raise RuntimeError(
        f"AUC矩阵细胞数为{auc.shape[0]}，"
        f"但merged_as有{merged_as.n_cells}个细胞。"
    )
if len(regulons) == 0:
    raise RuntimeError(
        "未得到regulon。请检查TF匹配、Hotspot FDR和motif pruning结果。"
    )
gc.collect()
import os
import numpy as np
import stereo as st
CLEAN_H5MS = "/data/work/GRN/AS_multi_GRN/AS_MSdata_GRN_clean.h5ms"
def remove_stale_raw(data, label):
    adata = data._ann_data
    print(f"\n[{label}]")
    print("current shape:", adata.shape)
    if adata.raw is not None:
        print("raw shape before removal:", adata.raw.shape)
        adata.raw = None
    else:
        print("raw already None")
    if hasattr(data, "_tl") and hasattr(data._tl, "_raw"):
        data._tl._raw = None
    if adata.raw is not None:
        raise RuntimeError(f"{label}: raw清除失败。")
    if data.exp_matrix.shape != (data.n_cells, data.n_genes):
        raise RuntimeError(
            f"{label}: exp_matrix维度与n_cells/n_genes不一致。"
        )
    for key, value in adata.obsm.items():
        if value.shape[0] != data.n_cells:
            raise RuntimeError(
                f"{label}: obsm[{key}]有{value.shape[0]}行，"
                f"但对象有{data.n_cells}个细胞。"
            )
    print("raw removed:", adata.raw is None)
    print("current obsm:", list(adata.obsm.keys()))
for sample_name in as_ms_data.names:
    remove_stale_raw(
        as_ms_data[sample_name],
        sample_name
    )
if as_ms_data.merged_data is not None:
    remove_stale_raw(
        as_ms_data.merged_data,
        "merged_data"
    )
total_cells = 0
for sample_name in as_ms_data.names:
    data = as_ms_data[sample_name]
    total_cells += data.n_cells
    if data.raw is not None:
        raise RuntimeError(
            f"{sample_name}: StereoPy raw仍然存在。"
        )
    print(
        sample_name,
        data.n_cells,
        data.n_genes,
        data.exp_matrix.shape
    )
if as_ms_data.merged_data.n_cells != total_cells:
    raise RuntimeError(
        "merged_data细胞数与10个样本细胞数之和不一致。"
    )
print("\nTotal cells:", total_cells)
print(
    "Merged shape:",
    as_ms_data.merged_data.exp_matrix.shape
)
st.io.write_h5ms(
    as_ms_data,
    output=CLEAN_H5MS
)
print("\nSaved clean H5MS:")
print(CLEAN_H5MS)
