# 模块 A · 读入 10X Visium（低精度，多样本）   [Python]
# 对应课程：第3课
# 运行：在项目目录下 `python load_visium.py`，或在 Jupyter / VS Code 里逐个代码块运行

# %%
# 步骤：导入包
# 输入：无
# 输出：无
from pathlib import Path
import pandas as pd
import scanpy as sc
import yaml

# %%
# 步骤：读取配置和样本表
# 输入：config.yaml；samples.tsv（sample / condition / path，path 指向 Space Ranger 的 outs 目录）
# 输出：cfg（字典）、samples（表格）、out_dir
cfg = yaml.safe_load(open("config.yaml"))
out_dir = Path(cfg["out_dir"])
out_dir.mkdir(parents=True, exist_ok=True)
samples = pd.read_csv(cfg["samples"], sep="\t")
samples

# %%
# 步骤：逐个样本读入，记录分组，并把坐标换算成 µm
# 输入：samples
# 输出：adatas（样本名 -> AnnData），每个样本都有 obsm['spatial']（像素）和 obsm['spatial_um']（µm）
adatas = {}
for row in samples.itertuples():
    adata = sc.read_visium(row.path, library_id=row.sample)
    adata.var_names_make_unique()
    adata.obs["condition"] = row.condition
    scalefactors = adata.uns["spatial"][row.sample]["scalefactors"]
    um_per_px = 55 / scalefactors["spot_diameter_fullres"]    # Visium spot 直径为 55µm
    adata.obsm["spatial_um"] = adata.obsm["spatial"] * um_per_px
    adatas[row.sample] = adata

# %%
# 步骤：合并所有样本，保存原始 counts
# 输入：adatas
# 输出：adata（spot x 基因），obs['sample']；layers['counts']
adata = sc.concat(adatas, label="sample", index_unique="_", uns_merge="unique")
adata.layers["counts"] = adata.X.copy()
adata

# %%
# 步骤：保存
# 输入：adata
# 输出：results/01_raw.h5ad
adata.write_h5ad(out_dir / "01_raw.h5ad")
