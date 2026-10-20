# 模块 B · 空间蛋白组（CODEX / IMC / mIHC）的转换、标准化、聚类   [Python]
# 对应课程：第26课
# 说明：蛋白强度不能用 RNA 的 normalize_total + log1p；这里用 arcsinh 转换后再按 marker 标准化。marker 只有几十个，所以不做 PCA。

# %%
# 步骤：导入包
# 输入：无
# 输出：无
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import squidpy as sq
import yaml

# %%
# 步骤：读取配置，设定 arcsinh 系数
# 输入：config.yaml
# 输出：cfg、out_dir、cofactor
cfg = yaml.safe_load(open("config.yaml"))
out_dir = Path(cfg["out_dir"])
cofactor = 5          # IMC / CyTOF 常用 5；CODEX 荧光强度大，常用 150 左右，先看强度分布再定

# %%
# 步骤：读入数据，arcsinh 转换
# 输入：results/01_raw.h5ad（layers['raw'] 为原始强度）
# 输出：adata（X 为 arcsinh 后的强度）
adata = sc.read_h5ad(out_dir / "01_raw.h5ad")
adata.X = np.arcsinh(adata.layers["raw"] / cofactor)

# %%
# 步骤：去掉几乎没有变化的 marker（空通道），否则下一步标准化会得到 NaN
# 输入：adata
# 输出：adata（marker 更少），打印被去掉的 marker
dead = adata.var_names[adata.X.std(axis=0) < 1e-3]
adata = adata[:, ~adata.var_names.isin(dead)].copy()
list(dead)

# %%
# 步骤：按 marker 做 z-score（截断在 ±5）
# 输入：adata
# 输出：adata（X 为标准化后的强度）
sc.pp.scale(adata, max_value=5)
adata

# %%
# 步骤：直接在 marker 空间建 KNN 图、UMAP、Leiden 聚类
# 输入：adata
# 输出：obs['leiden']、obsm['X_umap']
sc.pp.neighbors(adata, use_rep="X", n_neighbors=15, random_state=cfg["seed"])
sc.tl.umap(adata, random_state=cfg["seed"])
sc.tl.leiden(adata, resolution=cfg["leiden_resolution"], flavor="igraph", n_iterations=2, random_state=cfg["seed"])
adata.obs["leiden"].value_counts()

# %%
# 步骤：画每个聚类的平均 marker 强度，用于人工注释
# 输入：adata
# 输出：图 results/marker_matrixplot.png
sc.pl.matrixplot(adata, adata.var_names, groupby="leiden", cmap="RdBu_r", vcenter=0, show=False)
plt.savefig(out_dir / "marker_matrixplot.png", bbox_inches="tight")

# %%
# 步骤：画空间分布
# 输入：adata
# 输出：图 results/spatial_leiden.png
sq.pl.spatial_scatter(adata, color="leiden", library_key="sample", spatial_key="spatial_um", shape=None)
plt.savefig(out_dir / "spatial_leiden.png", bbox_inches="tight")

# %%
# 步骤：保存
# 输入：adata
# 输出：results/02_clustered.h5ad
adata.write_h5ad(out_dir / "02_clustered.h5ad")
