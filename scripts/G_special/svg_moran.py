# 模块 G · 空间可变基因（SVG）：Moran's I   [Python]
# 对应课程：第7课（SpaGene 找空间模式基因）的 Python 替代
# 适用：所有平台（CODEX 上就是"空间可变蛋白"）
# 多样本：逐样本计算，再汇总。Moran's I 衡量"邻居之间是否相似"；如果多个样本合在一起算，
#         样本之间整体表达水平的差异（批次）也会让同一样本内的邻居显得相似，把批次基因误判成空间可变基因。
#         所以每个样本单独算，再按"在几个样本里显著"和平均 I 排序。

# 说明：squidpy 的置换检验都加了 show_progress_bar=False。进度条会启动一个 multiprocessing 子进程，
#      在 macOS / Windows 上用 `python 脚本.py` 运行时，子进程会把整个脚本重新执行一遍，然后报错 EOFError。

# %%
# 步骤：导入包
# 输入：无
# 输出：无
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import squidpy as sq
import yaml

# %%
# 步骤：读取配置和数据，确定要检验的基因
# 输入：config.yaml、results/05_niche.h5ad
# 输出：adata、genes（高变基因；靶向 panel 和蛋白组为全部基因 / marker）、n_top
cfg = yaml.safe_load(open("config.yaml"))
out_dir = Path(cfg["out_dir"])
adata = sc.read_h5ad(out_dir / "05_niche.h5ad")
genes = adata.var_names[adata.var["highly_variable"]] if "highly_variable" in adata.var else adata.var_names
n_top = 6

# %%
# 步骤：每个样本单独建邻居图（6 个最近邻），计算 Moran's I 和置换 p 值
# 输入：adata、genes
# 输出：per_sample（样本名 -> Moran 结果表）；每个样本一张表 results/svg_moran_<样本>.csv
per_sample = {}
for sample in adata.obs["sample"].unique():
    sub = adata[adata.obs["sample"] == sample].copy()
    sq.gr.spatial_neighbors(sub, spatial_key="spatial_um", coord_type="generic", n_neighs=6)
    sq.gr.spatial_autocorr(sub, mode="moran", genes=genes, n_perms=100, n_jobs=cfg["threads"], seed=cfg["seed"], show_progress_bar=False)
    per_sample[sample] = sub.uns["moranI"]
    per_sample[sample].to_csv(out_dir / f"svg_moran_{sample}.csv")

# %%
# 步骤：汇总各样本：平均 I，以及在几个样本里显著（FDR < 0.05）
# 输入：per_sample
# 输出：summary（按显著样本数、平均 I 排序）；表 results/svg_moran.csv
moran_i = pd.DataFrame({s: t["I"] for s, t in per_sample.items()})
significant = pd.DataFrame({s: t["pval_norm_fdr_bh"] < 0.05 for s, t in per_sample.items()})
summary = pd.DataFrame({"mean_I": moran_i.mean(axis=1), "n_samples_significant": significant.sum(axis=1), "n_samples": len(per_sample)})
summary = summary.sort_values(["n_samples_significant", "mean_I"], ascending=False)
summary.to_csv(out_dir / "svg_moran.csv")
summary.head(10)

# %%
# 步骤：画汇总后最空间化的几个基因
# 输入：summary、adata
# 输出：图 results/svg_top.png
sq.pl.spatial_scatter(adata, color=summary.index[:n_top].tolist(), library_key="sample", spatial_key="spatial_um", shape=None)
plt.savefig(out_dir / "svg_top.png", bbox_inches="tight")
