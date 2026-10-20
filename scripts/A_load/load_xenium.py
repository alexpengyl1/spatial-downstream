# 模块 A · 读入 10X Xenium（单细胞级，靶向 panel）   [Python]
# 对应课程：第19课

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
# 输入：config.yaml；samples.tsv（path 指向 Xenium 输出目录）
# 输出：cfg、samples、out_dir
cfg = yaml.safe_load(open("config.yaml"))
out_dir = Path(cfg["out_dir"])
out_dir.mkdir(parents=True, exist_ok=True)
samples = pd.read_csv(cfg["samples"], sep="\t")
samples

# %%
# 步骤：逐个样本读入细胞 x 基因矩阵和细胞中心坐标（Xenium 坐标本身就是 µm）
# 输入：cell_feature_matrix.h5、cells.parquet
# 输出：adatas（样本名 -> AnnData），obsm['spatial'] = obsm['spatial_um']
adatas = {}
for row in samples.itertuples():
    adata = sc.read_10x_h5(Path(row.path) / "cell_feature_matrix.h5")    # 只保留 Gene Expression，去掉对照探针
    cells = pd.read_parquet(Path(row.path) / "cells.parquet").set_index("cell_id")
    adata.obs = cells.loc[adata.obs_names]
    adata.obs["condition"] = row.condition
    adata.obsm["spatial_um"] = adata.obs[["x_centroid", "y_centroid"]].to_numpy()
    adata.obsm["spatial"] = adata.obsm["spatial_um"]
    adatas[row.sample] = adata

# %%
# 步骤：合并样本，保存原始 counts
# 输入：adatas
# 输出：adata（细胞 x 基因）；layers['counts']
adata = sc.concat(adatas, label="sample", index_unique="_")
adata.layers["counts"] = adata.X.copy()
adata

# %%
# 步骤：保存
# 输入：adata
# 输出：results/01_raw.h5ad
adata.write_h5ad(out_dir / "01_raw.h5ad")
