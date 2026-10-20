# 模块 G · 按微环境给某一类细胞分亚型   [Python]
# 对应课程：第24课（附的文献是 CAF 的空间亚型）
# 适用：单细胞级平台（HD 细胞分割、Xenium、CosMx、CODEX）
# 思路：只取目标细胞（比如成纤维细胞），用它们的邻域组成（niche_composition.py 算好的 obsm['nhood_composition']）做 NMF，
#      再聚类：同一类细胞，周围环境不同，就分成不同的"空间亚型"；最后找各亚型的差异基因，看是否有功能差别。

# %%
# 步骤：导入包
# 输入：无
# 输出：无
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
from sklearn.decomposition import NMF
import squidpy as sq
import yaml

# %%
# 步骤：读取配置，指定目标细胞类型
# 输入：config.yaml
# 输出：cfg、out_dir、target_cell_type、n_components、subtype_resolution
cfg = yaml.safe_load(open("config.yaml"))
out_dir = Path(cfg["out_dir"])
target_cell_type = "Fibroblasts"
n_components = 5
subtype_resolution = 0.3

# %%
# 步骤：读入数据，取目标细胞
# 输入：results/05_niche.h5ad
# 输出：target（只含目标细胞）
adata = sc.read_h5ad(out_dir / "05_niche.h5ad")
target = adata[adata.obs["cell_type"] == target_cell_type].copy()
target

# %%
# 步骤：对目标细胞的邻域组成做 NMF
# 输入：target.obsm['nhood_composition']
# 输出：target.obsm['X_nhood_nmf']；loadings（成分 x 邻居细胞类型）
composition = pd.DataFrame(target.obsm["nhood_composition"], index=target.obs_names)
nmf = NMF(n_components=n_components, init="nndsvda", max_iter=500, random_state=cfg["seed"])
target.obsm["X_nhood_nmf"] = nmf.fit_transform(composition)
loadings = pd.DataFrame(nmf.components_, columns=composition.columns).round(2)
loadings

# %%
# 步骤：在 NMF 空间聚类，得到空间亚型
# 输入：target
# 输出：target.obs['spatial_subtype']
sc.pp.neighbors(target, use_rep="X_nhood_nmf", random_state=cfg["seed"])
sc.tl.leiden(target, resolution=subtype_resolution, key_added="spatial_subtype", flavor="igraph", n_iterations=2, random_state=cfg["seed"])
target.obs["spatial_subtype"].value_counts()

# %%
# 步骤：各亚型周围的细胞组成，以及在各分组中的占比
# 输入：target
# 输出：图 results/subtype_nhood.png；表 results/subtype_by_condition.csv
subtype_means = composition.groupby(target.obs["spatial_subtype"], observed=True).mean()
subtype_means.plot(kind="bar", stacked=True)
plt.savefig(out_dir / "subtype_nhood.png", bbox_inches="tight")
pd.crosstab(target.obs[cfg["condition_key"]], target.obs["spatial_subtype"], normalize="index").to_csv(out_dir / "subtype_by_condition.csv")

# %%
# 步骤：亚型之间的差异基因
# 输入：target（X 为 log 标准化表达）
# 输出：表 results/subtype_markers.csv
sc.tl.rank_genes_groups(target, groupby="spatial_subtype", method="wilcoxon", corr_method="benjamini-hochberg")
markers = sc.get.rank_genes_groups_df(target, group=None)
markers.to_csv(out_dir / "subtype_markers.csv", index=False)
markers.groupby("group", observed=True).head(5)

# %%
# 步骤：画亚型的空间分布
# 输入：target
# 输出：图 results/spatial_subtype.png
sq.pl.spatial_scatter(target, color="spatial_subtype", library_key="sample", spatial_key="spatial_um", shape=None)
plt.savefig(out_dir / "spatial_subtype.png", bbox_inches="tight")
