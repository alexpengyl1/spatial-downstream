# 模块 A · 读入 10X Visium HD 的方形 bin（高精度）   [Python]
# 对应课程：第11课
# 说明：读 binned_outputs/square_016um（或 008um）；如果用 Space Ranger 的细胞分割结果（segmented_outputs），
#      改读其中的 filtered_feature_cell_matrix.h5，坐标用 cell_segmentations.geojson 中多边形的中心点（见 references/platforms.md）。

# %%
# 步骤：导入包
# 输入：无
# 输出：无
from pathlib import Path
import json
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import yaml

# %%
# 步骤：读取配置和样本表
# 输入：config.yaml；samples.tsv（path 指向 Space Ranger 的 outs 目录）
# 输出：cfg、samples、out_dir、bin_name
cfg = yaml.safe_load(open("config.yaml"))
out_dir = Path(cfg["out_dir"])
out_dir.mkdir(parents=True, exist_ok=True)
samples = pd.read_csv(cfg["samples"], sep="\t")
bin_name = f"square_{cfg['hd_bin_um']:03d}um"     # 例如 square_016um
samples

# %%
# 步骤：逐个样本读入 bin x 基因的表达矩阵
# 输入：outs/binned_outputs/square_XXXum/filtered_feature_bc_matrix.h5
# 输出：adatas（样本名 -> AnnData），obs['condition']
adatas = {}
for row in samples.itertuples():
    bin_dir = Path(row.path) / "binned_outputs" / bin_name
    adatas[row.sample] = sc.read_10x_h5(bin_dir / "filtered_feature_bc_matrix.h5")
    adatas[row.sample].var_names_make_unique()
    adatas[row.sample].obs["condition"] = row.condition

# %%
# 步骤：给每个样本加上 bin 坐标（像素和 µm）
# 输入：spatial/tissue_positions.parquet、spatial/scalefactors_json.json
# 输出：adatas 中每个样本的 obsm['spatial']（像素）、obsm['spatial_um']（µm）
for row in samples.itertuples():
    spatial_dir = Path(row.path) / "binned_outputs" / bin_name / "spatial"
    positions = pd.read_parquet(spatial_dir / "tissue_positions.parquet").set_index("barcode")
    scalefactors = json.load(open(spatial_dir / "scalefactors_json.json"))
    adata = adatas[row.sample]
    adata.obsm["spatial"] = positions.loc[adata.obs_names, ["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy()
    adata.obsm["spatial_um"] = adata.obsm["spatial"] * scalefactors["microns_per_pixel"]

# %%
# 步骤：给每个样本加上 HE 图像和缩放系数，供画图时叠加
# 输入：spatial/tissue_hires_image.png、spatial/scalefactors_json.json
# 输出：adatas 中每个样本的 uns['spatial']
for row in samples.itertuples():
    spatial_dir = Path(row.path) / "binned_outputs" / bin_name / "spatial"
    hires = plt.imread(spatial_dir / "tissue_hires_image.png")
    scalefactors = json.load(open(spatial_dir / "scalefactors_json.json"))
    adatas[row.sample].uns["spatial"] = {row.sample: {"images": {"hires": hires}, "scalefactors": scalefactors}}

# %%
# 步骤：合并样本，保存原始 counts
# 输入：adatas
# 输出：adata（bin x 基因），obs['sample']；layers['counts']
adata = sc.concat(adatas, label="sample", index_unique="_", uns_merge="unique")
adata.layers["counts"] = adata.X.copy()
adata

# %%
# 步骤：保存
# 输入：adata
# 输出：results/01_raw.h5ad
adata.write_h5ad(out_dir / "01_raw.h5ad")
