import os
from natsort import natsorted
import stereo as st
from stereo.core.ms_data import MSData
from stereo.core.ms_pipeline import slice_generator
import warnings
warnings.filterwarnings('ignore')
import os
import scanpy as sc
import stereo as st
from natsort import natsorted
from stereo.core.ms_data import MSData
BASE_DIR = "/data/work/MS_subtype_NEW/H5ad_file/Region_noGM"
NOGM_DIR = "/data/work/MS_subtype_NEW/H5ad_file/Region_noGM"
NOGM_SAMPLES = {
    "Dr5_A5B5",
    "Dr6_A5B5",
    "Dr8_E5F6",
}
base_files = natsorted([
    os.path.join(BASE_DIR, f)
    for f in os.listdir(BASE_DIR)
    if f.endswith(".h5ad")
])
files = []
for base_path in base_files:
    sample_base = os.path.basename(base_path)
    sample_name = (
        sample_base
        .replace("filtered_", "")
        .replace("_subtype.h5ad", "")
        .replace(".h5ad", "")
    )
    if sample_name in NOGM_SAMPLES:
        sample_path = os.path.join(NOGM_DIR, sample_base)
    else:
        sample_path = base_path
    if not os.path.exists(sample_path):
        raise FileNotFoundError(sample_path)
    files.append(sample_path)
print("Input files:")
for f in files:
    print(f)
ms_data = MSData(_relationship="other", _var_type="intersect")
for sample_path in files:
    sample_base = os.path.basename(sample_path)
    sample_name = (
        sample_base
        .replace("filtered_", "")
        .replace("_subtype.h5ad", "")
        .replace(".h5ad", "")
    )
    adata = sc.read_h5ad(sample_path)
    if adata.raw is None:
        raise ValueError(f"raw is None: {sample_name}")
    if adata.raw.X.shape[1] != adata.n_vars:
        raise ValueError(
            f"raw.X genes != adata.n_vars in {sample_name}: "
            f"raw={adata.raw.X.shape[1]}, adata={adata.n_vars}"
        )
    if "X_spatial" not in adata.obsm:
        raise ValueError(f"X_spatial not found in adata.obsm: {sample_name}")
    if "subtype" not in adata.obs.columns:
        raise ValueError(f"subtype not found in adata.obs: {sample_name}")
    gene_names = adata.var_names.astype(str).tolist()
    adata.X = adata.raw.X.copy()
    adata.var_names = gene_names
    adata.var["gene_name"] = gene_names
    adata.var_names_make_unique()
    stereo_data = st.io.read_h5ad(
        anndata=adata,
        spatial_key="X_spatial",
        bin_type="cell_bins"
    )
    ms_data += stereo_data
ms_data.names = [
    os.path.basename(f)
    .replace("filtered_", "")
    .replace("_subtype.h5ad", "")
    .replace(".h5ad", "")
    for f in files
]
ms_data.integrate()
ms_data
from stereo.core.ms_pipeline import slice_generator
ms_data.to_integrate(
    res_key='subtype',
    scope=slice_generator[:],
    _from=slice_generator[:],
    type='obs',
    item=['subtype'] * ms_data.num_slice
)
print(ms_data.merged_data.tl.result['subtype'])
ms_data.tl.gen_ccc_micro_envs(
            cluster_res_key='subtype',
            n_boot=20,
            boot_prop=0.8,
            dimension=2,
            fill_rare=True,
            min_num=30,
            show_dividing_by_thresholds=True,
            method='split',
            threshold=None,
            scope=slice_generator[:],
            mode='integrate',
            output_path='/data/work/Stereopy_CCC_NEW/WM_region_subtype/',
            res_key='ccc_micro_envs',
            )
ms_data.tl.gen_ccc_micro_envs(
            scope=slice_generator[:],
            mode='integrate',
            method='split',
            threshold=0.186,
            output_path='/data/work/Stereopy_CCC_NEW/WM_region_subtype/micro_envs_2026-07-20_19-38-25',
            res_key='ccc_micro_envs'
            )
print(ms_data.merged_data.tl.result['ccc_micro_envs']['micro_envs'])
ms_data.tl.normalize_total(target_sum=None)
ms_data.tl.cell_cell_communication(
        scope=slice_generator[:],
        mode='integrate',
        analysis_type='statistical',
        cluster_res_key='subtype',
        species='HUMAN', 
        database='liana',  
        micro_envs='ccc_micro_envs',
        res_key='cell_cell_communication'
        )
data = ms_data.merged_data.adata
print("uns keys:", list(data.uns.keys()))
print("tl result keys:", ms_data.merged_data.tl.result.keys())
import os
import pandas as pd
import pickle
out = "/data/work/Stereopy_CCC_NEW/WM_region_subtype/CCC_export/"
os.makedirs(out, exist_ok=True)
ccc = ms_data.merged_data.adata.uns["cell_cell_communication"]
if hasattr(ccc, "to_csv"):
    ccc.to_csv(out + "cell_cell_communication.csv")
else:
    with open(out + "cell_cell_communication.pkl", "wb") as f:
        pickle.dump(ccc, f)
micro = ms_data.merged_data.adata.uns["ccc_micro_envs"]
if hasattr(micro, "to_csv"):
    micro.to_csv(out + "microenv.csv")
else:
    with open(out + "microenv.pkl", "wb") as f:
        pickle.dump(micro, f)
print("done")
import pandas as pd
import stereo as st
data = ms_data.merged_data
sample_names = [
'Dr1_A2_region_filtered', 'Dr1_F4_region_filtered', 'Dr2_D6_region_filtered', 
    'Dr3_F4_region_filtered', 'Dr4_F2_region_filtered', 'Dr5_A5B5_region_filtered', 'Dr6_A2_region_filtered', 
    'Dr6_A5B5_region_filtered', 'Dr8_E3_region_filtered', 'Dr8_E5F6_region_filtered'
]
if list(ms_data.names) != sample_names:
    raise ValueError(f"ms_data.names mismatch:\n{ms_data.names}")
data.adata.uns["sn"] = pd.DataFrame({
    "batch": [str(i) for i in range(len(sample_names))],
    "sn": sample_names
})
print(data.adata.uns["sn"])
OUT_H5MS = "/data/work/Stereopy_CCC_NEW/WM_region_subtype/MS—subtype_microenv_with_CCC.h5ms"
for i, sample_name in enumerate(sample_names):
    sd = ms_data[i]
    sad = sd.adata
    sad.uns["sn"] = pd.DataFrame({
        "batch": [str(i)],
        "sn": [sample_name]
    })
    print(i, sample_name, "sn type:", type(sad.uns["sn"]), sad.uns["sn"].shape)
    assert isinstance(sad.uns["sn"], pd.DataFrame)
st.io.write_h5ms(
    ms_data,
    output=OUT_H5MS
)
ms_data.merged_data.tl.result['cell_cell_communication']['significant_means']
significant_means=ms_data.merged_data.tl.result['cell_cell_communication']['significant_means']
output_path = '/data/work/Stereopy_CCC_NEW/WM_region_subtype/AS-subtype_CCC_significant_means.csv'
significant_means.to_csv(output_path, index=False)
