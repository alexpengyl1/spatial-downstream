# 模块 G · 空间 CNV（区分肿瘤 / 正常区域，看克隆的空间分布）   [Python]
# 对应课程：第2、8课（课程第2课用 R 的 inferCNV；infercnvpy 是它的 Python 实现，结果可比，速度快得多）
# 适用：全转录组平台（Visium、HD、Stereo-seq）。靶向 panel（Xenium / CosMx / CODEX）基因太少，不要做。
# 前提：config 里填好 gtf（用来取基因的染色体位置）和 cnv_reference（正常参考细胞类型，比如免疫细胞）。
# HD：2/8µm bin 太稀疏，用 16µm bin 或细胞分割结果。

# %%
# 步骤：导入包
# 输入：无
# 输出：无
from pathlib import Path
import multiprocessing
import infercnvpy as cnv
import matplotlib.pyplot as plt
import scanpy as sc
import squidpy as sq
import yaml

# infercnvpy 一定会开进程池；macOS / Windows 默认用 spawn，用 `python 脚本.py` 运行时子进程会重新执行整个脚本并报错。
# 改用 fork 就没有这个问题（Windows 不支持 fork，请删掉这一行并在 Jupyter 里运行）
multiprocessing.set_start_method("fork", force=True)

# %%
# 步骤：读取配置
# 输入：config.yaml
# 输出：cfg、out_dir、window_size
cfg = yaml.safe_load(open("config.yaml"))
out_dir = Path(cfg["out_dir"])
window_size = 250          # 平滑窗口包含的基因数；infercnvpy 默认 100，空间数据噪声大，用 250

# %%
# 步骤：读入数据，从 GTF 取基因的染色体位置
# 输入：results/05_niche.h5ad、gtf
# 输出：adata.var 里多了 chromosome / start / end
adata = sc.read_h5ad(out_dir / "05_niche.h5ad")
cnv.io.genomic_position_from_gtf(cfg["gtf"], adata)
adata.var["chromosome"].value_counts().head()

# %%
# 步骤：以正常参考细胞为基线推断 CNV
# 输入：adata（X 为 log 标准化表达）、cnv_reference
# 输出：obsm['X_cnv']
cnv.tl.infercnv(adata, reference_key="cell_type", reference_cat=cfg["cnv_reference"], window_size=window_size)

# %%
# 步骤：在 CNV 矩阵上聚类，并计算每个 spot / 细胞的 CNV 得分
# 输入：obsm['X_cnv']
# 输出：obs['cnv_leiden']、obs['cnv_score']
cnv.tl.pca(adata)
cnv.pp.neighbors(adata)
cnv.tl.leiden(adata, resolution=0.3)
cnv.tl.cnv_score(adata)
adata.obs.groupby("cnv_leiden", observed=True)["cnv_score"].mean()

# %%
# 步骤：画按 CNV 聚类分组的染色体热图
# 输入：adata
# 输出：图 results/cnv_heatmap.png
adata.obsm["X_cnv"] = adata.obsm["X_cnv"].toarray()      # 稀疏矩阵在新版 pandas / scanpy 下画热图会报错
cnv.pl.chromosome_heatmap(adata, groupby="cnv_leiden", show=False)
plt.savefig(out_dir / "cnv_heatmap.png", bbox_inches="tight")

# %%
# 步骤：画 CNV 得分和 CNV 聚类的空间分布
# 输入：adata
# 输出：图 results/spatial_cnv.png
sq.pl.spatial_scatter(adata, color=["cnv_score", "cnv_leiden"], library_key="sample", spatial_key="spatial_um", shape=None)
plt.savefig(out_dir / "spatial_cnv.png", bbox_inches="tight")

# %%
# 步骤：保存 CNV 结果表
# 输入：adata
# 输出：表 results/cnv_per_spot.csv
adata.obs[["sample", "cell_type", "niche", "cnv_leiden", "cnv_score"]].to_csv(out_dir / "cnv_per_spot.csv")
