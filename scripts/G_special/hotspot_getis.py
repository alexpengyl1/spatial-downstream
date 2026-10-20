# 模块 G · 热点分析（Getis-Ord Gi*）和热点之间的距离   [Python]
# 对应课程：第9课（课程用 spottedpy，它底层也是 Getis-Ord；这里直接用 PySAL 的 esda，便于控制参数）
# 适用：任何"每个 spot / 细胞一个数值"的特征：某类细胞的比例、基因集得分、CNV 得分 ……
# 问题示例：成纤维细胞的热点区，离增殖肿瘤细胞的热点区有多远？

# %%
# 步骤：导入包
# 输入：无
# 输出：无
from pathlib import Path
from esda.getisord import G_Local
from libpysal.weights import KNN
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.spatial import cKDTree
import squidpy as sq
import yaml

# %%
# 步骤：读取配置，指定两个要比较的特征
# 输入：config.yaml
# 输出：cfg、out_dir、feature_a、feature_b、k
cfg = yaml.safe_load(open("config.yaml"))
out_dir = Path(cfg["out_dir"])
feature_a = "Fibroblasts"     # 细胞类型名（取自细胞组成矩阵）或 obs 里的数值列名
feature_b = "Epithelial"
k = 10                        # 每个点用几个最近邻定义"局部"

# %%
# 步骤：读入数据，把细胞组成的各列放进 obs，方便按名字取特征
# 输入：results/05_niche.h5ad
# 输出：adata
adata = sc.read_h5ad(out_dir / "05_niche.h5ad")
key = cfg["composition_key"] if cfg["resolution"] == "spot" else "nhood_composition"
composition = pd.DataFrame(adata.obsm[key], index=adata.obs_names)
adata.obs[composition.columns] = composition

# %%
# 步骤：每个样本单独计算 Gi*，z > 0 且 p < 0.05 为热点（hot），z < 0 且 p < 0.05 为冷点（cold）
# 输入：adata
# 输出：obs['<特征>_hotspot']（hot / cold / ns）
for feature in (feature_a, feature_b):
    adata.obs[f"{feature}_hotspot"] = "ns"
    for sample in adata.obs["sample"].unique():
        in_sample = (adata.obs["sample"] == sample).to_numpy()
        weights = KNN.from_array(adata.obsm["spatial_um"][in_sample], k=k)
        gi = G_Local(adata.obs.loc[in_sample, feature].to_numpy(), weights, star=True, permutations=999, seed=cfg["seed"])
        label = np.where(gi.p_sim < 0.05, np.where(gi.Zs > 0, "hot", "cold"), "ns")
        adata.obs.loc[in_sample, f"{feature}_hotspot"] = label
    adata.obs[f"{feature}_hotspot"] = adata.obs[f"{feature}_hotspot"].astype("category")
adata.obs[[f"{feature_a}_hotspot", f"{feature_b}_hotspot"]].value_counts()

# %%
# 步骤：画两个特征的热点分布
# 输入：adata
# 输出：图 results/hotspot.png
sq.pl.spatial_scatter(adata, color=[f"{feature_a}_hotspot", f"{feature_b}_hotspot"], library_key="sample", spatial_key="spatial_um", shape=None)
plt.savefig(out_dir / "hotspot.png", bbox_inches="tight")

# %%
# 步骤：计算 A 的热点和冷点到最近的 B 热点的距离（µm），每个样本单独算
# 输入：adata
# 输出：distances（表：样本、A 的热/冷、到 B 热点的距离）；表 results/hotspot_distance.csv
rows = []
for sample in adata.obs["sample"].unique():
    sub = adata[adata.obs["sample"] == sample]
    tree = cKDTree(sub.obsm["spatial_um"][(sub.obs[f"{feature_b}_hotspot"] == "hot").to_numpy()])
    for state in ("hot", "cold"):
        distance, _ = tree.query(sub.obsm["spatial_um"][(sub.obs[f"{feature_a}_hotspot"] == state).to_numpy()])
        rows.append(pd.DataFrame({"sample": sample, f"{feature_a}": state, f"distance_to_{feature_b}_hot": distance}))
distances = pd.concat(rows)
distances.to_csv(out_dir / "hotspot_distance.csv", index=False)
distances.groupby(["sample", feature_a])[f"distance_to_{feature_b}_hot"].median()

# %%
# 步骤：比较 A 的热点和冷点到 B 热点的距离
# 输入：distances
# 输出：图 results/hotspot_distance.png
distances.boxplot(column=f"distance_to_{feature_b}_hot", by=feature_a)
plt.savefig(out_dir / "hotspot_distance.png", bbox_inches="tight")
