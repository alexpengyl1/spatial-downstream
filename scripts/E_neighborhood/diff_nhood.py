# 模块 E · 两组之间的邻域富集差异   [Python]
# 对应课程：第16、23课
# 适用：单细胞级平台，两个分组（condition），每组至少 2 个样本（≥ 3 个更可靠）。
# 做法：以样本为重复单位：每个样本单独做邻域富集（z 分数），再对每一对细胞类型在两组之间做 Welch t 检验。
# 和课程代码的差别：
#   第16课用 co_occurrence 的第一个距离区间直接相减，没有统计检验；
#   第23课用 CellCharter 的 diff_nhood_enrichment，但 0.3.5 版在 n_jobs=1 时报参数错误；
#   n_jobs>1 时它会开 ProcessPoolExecutor，用 `python 脚本.py` 运行时在 macOS 上会卡住。

# %%
# 步骤：导入包
# 输入：无
# 输出：无
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import ttest_ind
import seaborn as sns
import squidpy as sq
import yaml

# %%
# 步骤：读取配置和 niche 结果，确定要比较的两组
# 输入：config.yaml、results/05_niche.h5ad
# 输出：adata、cluster_key、condition_key、group1 / group2；打印每组样本数
cfg = yaml.safe_load(open("config.yaml"))
out_dir = Path(cfg["out_dir"])
adata = sc.read_h5ad(out_dir / "05_niche.h5ad")
cluster_key = "cell_type"
condition_key = cfg["condition_key"]
adata.obs[condition_key] = adata.obs[condition_key].astype("category")
group1, group2 = adata.obs[condition_key].cat.categories[:2]      # 结果为 group2 - group1
adata.obs.groupby(condition_key, observed=True)["sample"].nunique()

# %%
# 步骤：每个样本单独建 Delaunay 邻居图、做邻域富集，收集 z 分数
# 输入：adata（obs['cell_type']、obsm['spatial_um']）
# 输出：zscores（列表，每个样本一个 细胞类型 x 细胞类型 的矩阵）
samples = adata.obs["sample"].unique()
zscores = []
for sample in samples:
    sub = adata[adata.obs["sample"] == sample].copy()
    sq.gr.spatial_neighbors(sub, spatial_key="spatial_um", coord_type="generic", delaunay=True)
    sq.gr.nhood_enrichment(sub, cluster_key=cluster_key, n_perms=1000, seed=cfg["seed"], show_progress_bar=False)
    zscores.append(sub.uns[f"{cluster_key}_nhood_enrichment"]["zscore"])

# %%
# 步骤：把 z 分数堆成三维数组，并记下每个样本属于哪一组
# 输入：zscores、samples
# 输出：zscores（样本 x 细胞类型 x 细胞类型）、conditions、cell_types
zscores = np.stack(zscores)
conditions = adata.obs.groupby("sample", observed=True)[condition_key].first().loc[samples].to_numpy()
cell_types = adata.obs[cluster_key].cat.categories
zscores.shape

# %%
# 步骤：两组之间比较每一对细胞类型的 z 分数（Welch t 检验）
# 输入：zscores、conditions
# 输出：diff（平均 z 分数之差，group2 - group1）、pvalues
in_group1, in_group2 = conditions == group1, conditions == group2
diff = pd.DataFrame(zscores[in_group2].mean(axis=0) - zscores[in_group1].mean(axis=0), index=cell_types, columns=cell_types)
pvalues = pd.DataFrame(ttest_ind(zscores[in_group2], zscores[in_group1], axis=0, equal_var=False).pvalue, index=cell_types, columns=cell_types)
diff.round(1)

# %%
# 步骤：保存差异和 p 值，画差异热图（* 表示 p < 0.05）
# 输入：diff、pvalues
# 输出：表 results/diff_nhood_zscore.csv、results/diff_nhood_pvalue.csv；图 results/diff_nhood.png
diff.to_csv(out_dir / "diff_nhood_zscore.csv")
pvalues.to_csv(out_dir / "diff_nhood_pvalue.csv")
sns.heatmap(diff, cmap="bwr", center=0, annot=np.where(pvalues < 0.05, "*", ""), fmt="", linewidths=0.5)
plt.title(f"nhood enrichment z: {group2} - {group1}")
plt.savefig(out_dir / "diff_nhood.png", bbox_inches="tight")
