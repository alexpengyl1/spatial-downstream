# 模块 D · 细胞 niche / RCN / CN：按"周围有哪些细胞"聚类   [Python]
# 对应课程：第5、14、21、27课
# 低精度平台（resolution: spot）：一个 spot 的细胞组成（obsm[composition_key]，来自 RCTD / cell2location / marker 打分）
#   本身就是一个小邻域，直接对它聚类（第5课）。
# 单细胞级平台（resolution: cell）：统计每个细胞半径 radius_um 内各类细胞的比例，再聚类（第14、21、27课的 RCN / CN）。
# 结果：obs['niche']；obsm['nhood_composition']（第24课的邻域亚型分析会用到）
# 多样本：联合。邻域组成按样本分别统计（邻居图不跨样本），KMeans 在所有样本上一起做，niche 编号跨样本可比；
#         另外按样本输出各 niche 的占比和细胞组成（第5课按病人取中位数再比较），供以样本为单位做统计。

# %%
# 步骤：导入包
# 输入：无
# 输出：无
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns
from sklearn.cluster import KMeans
import squidpy as sq
import yaml

# %%
# 步骤：读取配置；有空间域结果就读 04_domain，否则读 03_annotated
# 输入：config.yaml、results/04_domain.h5ad 或 03_annotated.h5ad
# 输出：cfg、out_dir、adata
cfg = yaml.safe_load(open("config.yaml"))
out_dir = Path(cfg["out_dir"])
input_file = out_dir / "04_domain.h5ad" if (out_dir / "04_domain.h5ad").exists() else out_dir / "03_annotated.h5ad"
adata = sc.read_h5ad(input_file)
adata

# %%
# 步骤：得到每个 spot/细胞的邻域组成矩阵
# 输入：低精度：obsm[composition_key]；单细胞级：obs['cell_type'] + 半径邻居图
# 输出：composition（DataFrame，行 = spot/细胞，列 = 细胞类型，每行和为 1）
if cfg["resolution"] == "spot":
    composition = pd.DataFrame(adata.obsm[cfg["composition_key"]], index=adata.obs_names).clip(lower=0)
else:
    sq.gr.spatial_neighbors(adata, library_key="sample", spatial_key="spatial_um", coord_type="generic", radius=cfg["radius_um"], set_diag=True)
    onehot = pd.get_dummies(adata.obs["cell_type"]).astype(float)
    composition = pd.DataFrame(adata.obsp["spatial_connectivities"] @ onehot.to_numpy(), index=adata.obs_names, columns=onehot.columns)
composition = composition.div(composition.sum(axis=1), axis=0).fillna(0)
adata.obsm["nhood_composition"] = composition
composition.head()

# %%
# 步骤：用 KMeans 把邻域组成聚成 n_niches 类
# 输入：composition
# 输出：obs['niche']（niche_1、niche_2 ...）
kmeans = KMeans(n_clusters=cfg["n_niches"], random_state=cfg["seed"], n_init=10)
adata.obs["niche"] = pd.Categorical([f"niche_{i + 1}" for i in kmeans.fit_predict(composition)])
adata.obs["niche"].value_counts()

# %%
# 步骤：每个 niche 的平均细胞组成（热图，按列标准化），用来给 niche 命名
# 输入：adata、composition
# 输出：图 results/niche_composition_heatmap.png
niche_means = composition.groupby(adata.obs["niche"], observed=True).mean()
sns.clustermap(niche_means, standard_scale=1, cmap="Reds", figsize=(8, 5))
plt.savefig(out_dir / "niche_composition_heatmap.png", bbox_inches="tight")
niche_means.round(2)

# %%
# 步骤：每个样本里各 niche 的占比（比较不同分组用）
# 输入：adata
# 输出：表 results/niche_by_sample.csv；图 results/niche_by_sample.png
niche_by_sample = pd.crosstab(adata.obs["sample"], adata.obs["niche"], normalize="index")
niche_by_sample.to_csv(out_dir / "niche_by_sample.csv")
niche_by_sample.plot(kind="bar", stacked=True)
plt.savefig(out_dir / "niche_by_sample.png", bbox_inches="tight")

# %%
# 步骤：每个样本、每个 niche 的细胞组成中位数（以样本为单位比较 niche 的特征细胞时用）
# 输入：composition、adata.obs['sample' / 'niche']
# 输出：表 results/niche_composition_by_sample.csv（行 = 样本 × niche，列 = 细胞类型）
composition_by_sample = composition.groupby([adata.obs["sample"], adata.obs["niche"]], observed=True).median()
composition_by_sample.to_csv(out_dir / "niche_composition_by_sample.csv")
composition_by_sample.head()

# %%
# 步骤：画 niche 的空间分布
# 输入：adata
# 输出：图 results/spatial_niche.png
sq.pl.spatial_scatter(adata, color="niche", library_key="sample", spatial_key="spatial_um", shape=None)
plt.savefig(out_dir / "spatial_niche.png", bbox_inches="tight")

# %%
# 步骤：保存（后面的邻域、通讯、个性化分析都从这个文件开始）
# 输入：adata
# 输出：results/05_niche.h5ad
adata.write_h5ad(out_dir / "05_niche.h5ad")
