# 模块 E · 细胞类型的邻域富集和共现   [Python]
# 对应课程：第15、19、22课
# 适用：单细胞级平台（HD、Xenium、CosMx、CODEX）。低精度 Visium 改用 misty_colocalization.py（spot 里本来就是混合的）。
# 两个指标：
#   nhood_enrichment  直接相邻（Delaunay 连边）的两类细胞，比随机打乱标签多还是少（置换检验得到 z 分数）
#   co_occurrence     随距离变化的共现概率比，能看出两类细胞在多远的范围内聚在一起
# 多样本：逐样本计算，再汇总。置换检验打乱的是细胞类型标签，如果把所有样本合在一起打乱，
#         各样本细胞组成的差异会混进零分布；所以每个样本单独做（和课程第19课按 FOV 循环一致），再取各样本 z 分数的平均。
# 说明：squidpy 的置换检验都加了 show_progress_bar=False。进度条会启动一个 multiprocessing 子进程，
#      在 macOS / Windows 上用 `python 脚本.py` 运行时，子进程会把整个脚本重新执行一遍，然后报错 EOFError。

# %%
# 步骤：导入包
# 输入：无
# 输出：无
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import squidpy as sq
import yaml

# %%
# 步骤：读取配置和 niche 结果
# 输入：config.yaml、results/05_niche.h5ad
# 输出：cfg、out_dir、adata、cluster_key、cell_types
cfg = yaml.safe_load(open("config.yaml"))
out_dir = Path(cfg["out_dir"])
adata = sc.read_h5ad(out_dir / "05_niche.h5ad")
cluster_key = "cell_type"
cell_types = adata.obs[cluster_key].cat.categories

# %%
# 步骤：每个样本单独建 Delaunay 邻居图、做邻域富集（置换检验）
# 输入：adata（obs['cell_type']、obsm['spatial_um']）
# 输出：zscores（样本名 -> 细胞类型 x 细胞类型的 z 分数表）；每个样本一张表 results/nhood_enrichment_<样本>.csv
zscores = {}
for sample in adata.obs["sample"].unique():
    sub = adata[adata.obs["sample"] == sample].copy()
    sq.gr.spatial_neighbors(sub, spatial_key="spatial_um", coord_type="generic", delaunay=True)
    sq.gr.nhood_enrichment(sub, cluster_key=cluster_key, n_perms=1000, seed=cfg["seed"], show_progress_bar=False)
    zscores[sample] = pd.DataFrame(sub.uns[f"{cluster_key}_nhood_enrichment"]["zscore"], index=cell_types, columns=cell_types)
    zscores[sample].to_csv(out_dir / f"nhood_enrichment_{sample}.csv")

# %%
# 步骤：各样本 z 分数取平均（某样本里没有的细胞类型记为 NaN，不参与平均），并记录有多少个样本方向一致
# 输入：zscores
# 输出：mean_z、n_positive（z > 0 的样本数）；表 results/nhood_enrichment_mean.csv
stacked = np.stack([z.to_numpy() for z in zscores.values()])
mean_z = pd.DataFrame(np.nanmean(stacked, axis=0), index=cell_types, columns=cell_types)
n_positive = pd.DataFrame((stacked > 0).sum(axis=0), index=cell_types, columns=cell_types)
mean_z.to_csv(out_dir / "nhood_enrichment_mean.csv")
mean_z.round(1)

# %%
# 步骤：画平均 z 分数热图（对角线是同类细胞聚集，数值通常很大，设为空白，免得把其他格子的颜色压扁）；格子里写"z > 0 的样本数 / 样本总数"
# 输入：mean_z、n_positive
# 输出：图 results/nhood_enrichment_mean.png
plot_z = mean_z.copy()
np.fill_diagonal(plot_z.values, np.nan)
sns.heatmap(plot_z, cmap="bwr", center=0, annot=n_positive.astype(str) + f"/{len(zscores)}", fmt="", linewidths=0.5)
plt.title("mean nhood enrichment z (cells: samples with z > 0)")
plt.savefig(out_dir / "nhood_enrichment_mean.png", bbox_inches="tight")

# %%
# 步骤：随距离变化的共现（每个样本单独算；距离区间 0–300µm）
# 输入：adata
# 输出：每个样本一张图 results/co_occurrence_<样本>.png
for sample in adata.obs["sample"].unique():
    sub = adata[adata.obs["sample"] == sample].copy()
    sq.gr.co_occurrence(sub, cluster_key=cluster_key, spatial_key="spatial_um", interval=np.linspace(5, 300, 30), show_progress_bar=False)
    sq.pl.co_occurrence(sub, cluster_key=cluster_key, clusters=list(sub.obs[cluster_key].cat.categories[:4]))
    plt.savefig(out_dir / f"co_occurrence_{sample}.png", bbox_inches="tight")
