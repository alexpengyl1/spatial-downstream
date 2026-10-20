# 模块 A · 读入华大 Stereo-seq GEF   [Python，需要 stereopy 环境]
# 对应课程：第3、11课
# 说明：bin50 及以上按低精度（resolution: spot）处理；bin20 及以下或 cellbin 按高精度（resolution: cell）处理。
#      stereopy 依赖版本较老，建议单独建 stereo 环境，只用这个脚本把 GEF 转成 h5ad，后续步骤回到 spatial-py 环境。

# %%
# 步骤：导入包
# 输入：无
# 输出：无
from pathlib import Path
import pandas as pd
import scanpy as sc
import stereo as st
import yaml

# %%
# 步骤：读取配置和样本表
# 输入：config.yaml；samples.tsv（path 指向 .gef 文件）
# 输出：cfg、samples、out_dir
cfg = yaml.safe_load(open("config.yaml"))
out_dir = Path(cfg["out_dir"])
out_dir.mkdir(parents=True, exist_ok=True)
samples = pd.read_csv(cfg["samples"], sep="\t")
samples

# %%
# 步骤：逐个样本读入 GEF，转成 AnnData，坐标换算成 µm
# 输入：*.gef
# 输出：adatas（样本名 -> AnnData），obsm['spatial']、obsm['spatial_um']
adatas = {}
for row in samples.itertuples():
    stereo_data = st.io.read_gef(row.path, bin_size=cfg["stereo_bin_size"])    # cellbin 数据改为 bin_type="cell_bins"
    adata = st.io.stereo_to_anndata(stereo_data, flavor="scanpy")
    adata.var_names_make_unique()
    adata.obs["condition"] = row.condition
    adata.obsm["spatial_um"] = adata.obsm["spatial"] * cfg["stereo_um_per_unit"]
    adatas[row.sample] = adata

# %%
# 步骤：合并样本，保存原始 counts
# 输入：adatas
# 输出：adata；layers['counts']
adata = sc.concat(adatas, label="sample", index_unique="_")
adata.layers["counts"] = adata.X.copy()
adata

# %%
# 步骤：保存
# 输入：adata
# 输出：results/01_raw.h5ad
adata.write_h5ad(out_dir / "01_raw.h5ad")
