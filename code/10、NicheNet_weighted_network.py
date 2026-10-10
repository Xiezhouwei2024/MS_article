import os
from natsort import natsorted
import stereo as st
from stereo.core.ms_data import MSData
from stereo.core.ms_pipeline import slice_generator
import warnings
warnings.filterwarnings('ignore')
ms_data = st.io.read_h5ms(
    '/data/work/Stereopy_CCC_NEW/WM_region_subtype/MS—subtype_microenv_with_CCC.h5ms',
    use_raw=True,
    use_result=True
)
out_dir = "/data/work/GRN/NicheNe-V2-CCC-GRN-Sankey"
fig = ms_data.plt.ccc_sankey_plot(
    sender_cluster="MIMS-iron",
    receiver_cluster="AS_S6_Trans",
    weighted_network_path="/data/work/GRN/3D_GRN_human/NicheNet-V2_weighted_network_lr_sig_human_stereopy.tsv",
    regulons="/data/work/GRN/AS_multi_GRN/Final_use/AS_multi_hotspot_GRN_regulon_list.csv",
    pct_expressed=0.05,
    max_path_length=4,
    receptor_tf_paths_out_path=f"{out_dir}/MIMS-iron_AS-S6_receptor_TF_paths.csv",
    out_path=f"{out_dir}/MIMS-iron_AS-S6_CCC_GRN_sankey.png",
    width=1000,
    height=1000
)
