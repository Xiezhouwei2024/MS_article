import os
import gc
import warnings
from collections import defaultdict
import anndata as ad
import stereo as st
import numpy as np
import pandas as pd
from scipy import sparse
import matplotlib as mpl
warnings.filterwarnings("ignore")
ann_h5ad = "/Files/xiezhouwei/MS/Stereo_data/Bin50/MS_bin50_region_h5ad/Dr6_C02134A2_BIN50_region_modify.h5ad"
out_dir = "/data/work/Spatial_hotspot_bin50/Figures"
os.makedirs(out_dir, exist_ok=True)
sample_id = os.path.basename(ann_h5ad).replace("_region.h5ad", "")
print(sample_id)
hvg_key = "highly_variable_genes_4000"
hotspot_key = "spatial_hotspot"
adata = ad.read_h5ad(ann_h5ad)
if adata.raw is None:
    raise ValueError("adata.raw is missing; the original count matrix cannot be read.")
if "X_spatial" not in adata.obsm:
    raise ValueError("X_spatial is missing from adata.obsm.")
if adata.X.shape != adata.raw.X.shape:
    raise ValueError(
        f"X and raw.X have different shapes: X={adata.X.shape}, raw={adata.raw.X.shape}"
    )
if adata.var_names.has_duplicates:
    duplicated = adata.var_names[
        adata.var_names.duplicated()
    ].unique().tolist()
    raise ValueError(f"Duplicate gene names detected, for example: {duplicated[:20]}")
raw_test = adata.raw.X[:2000, :]
raw_values = (
    raw_test.data
    if sparse.issparse(raw_test)
    else np.asarray(raw_test)[np.asarray(raw_test) != 0]
)
integer_fraction = np.mean(
    np.isclose(raw_values, np.round(raw_values), atol=1e-8)
)
print("AnnData shape:", adata.shape)
print("raw exists:", adata.raw is not None)
print("spatial shape:", adata.obsm["X_spatial"].shape)
print("raw integer-like fraction:", integer_fraction)
print("raw nonzero min:", raw_values.min())
print("raw nonzero max:", raw_values.max())
if integer_fraction < 0.999:
    raise ValueError("raw.X is not an original integer count matrix; analysis stopped.")
data_bin50 = st.io.anndata_to_stereo(
    andata=adata,
    use_raw=True,
    spatial_key="X_spatial"
)
data_bin50.bin_type = "bins"
data_bin50.bin_size = 50
data_bin50.sn = sample_id
print("\nStereoExpData information")
print("shape:", data_bin50.exp_matrix.shape)
print("bin_type:", data_bin50.bin_type)
print("bin_size:", data_bin50.bin_size)
print("position shape:", data_bin50.position.shape)
print("first genes:", data_bin50.genes.gene_name[:20])
del adata, raw_test, raw_values
gc.collect()
data_bin50.tl.raw_checkpoint()
data_bin50.tl.normalize_total(
    target_sum=10000
)
data_bin50.tl.log1p()
data_bin50.tl.highly_variable_genes(
    min_mean=0.0125,
    max_mean=3,
    min_disp=0.5,
    n_top_genes=4000,
    res_key=hvg_key
)
data_bin50.plt.highly_variable_genes(
    res_key=hvg_key,
    out_path=os.path.join(
        out_dir,
        f"{sample_id}_HVG_4000.png"
    )
)
gc.collect()
data_bin50.tl.spatial_hotspot(
    use_highly_genes=True,
    use_raw=True,
    hvg_res_key=hvg_key,
    model="danb",
    n_neighbors=40,
    n_jobs=20,
    fdr_threshold=0.05,
    min_gene_threshold=10,
    res_key=hotspot_key
)
print("\nHotspot analysis completed")
print("result key:", hotspot_key)
mpl.rcParams["figure.dpi"] = 300
mpl.rcParams["savefig.dpi"] = 300
data_bin50.plt.hotspot_local_correlations(
    res_key="spatial_hotspot",
    out_path=os.path.join(
        out_dir,
        f"{sample_id}_heatmap.png"
    )
)
data_bin50.plt.hotspot_modules(
    res_key="spatial_hotspot",
    out_path=os.path.join(
        out_dir,
        f"{sample_id}_spatial.png"
    )
)
hotspot = data_bin50.tl.result["spatial_hotspot"]
module_gene_list = defaultdict(list)
for gene, module in hotspot.modules.items():
    if module == -1:
        continue
    module_gene_list[module].append(gene)
for module in sorted(module_gene_list):
    print(f"\nModule {module} ({len(module_gene_list[module])} genes)")
    print(module_gene_list[module])
out_file = os.path.join(
    out_dir,
    f"{sample_id}_Hotspot_module_genelist.csv"
)
rows = []
for module in sorted(module_gene_list):
    rows.append({
        "Module": module,
        "Gene_number": len(module_gene_list[module]),
        "Genes": ",".join(module_gene_list[module])
    })
pd.DataFrame(rows).to_csv(out_file, index=False)
print(f"\nSaved: {out_file}")
