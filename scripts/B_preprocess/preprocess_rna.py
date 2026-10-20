# 模块 B · RNA 平台的质控、标准化、去批次、聚类   [Python]
# 适用：Visium / Visium HD / Stereo-seq / Xenium / CosMx
# 对应课程：第1、3、11、19课
# 和课程代码的差别：去掉了 regress_out（大数据很慢，Harmony 已能处理批次）；靶向 panel 不选高变基因；保留全部基因，只在 PCA 时用高变基因。

# %%
# 步骤：导入包
# 输入：无
# 输出：无
from pathlib import Path
import matplotlib.pyplot as plt
import scanpy as sc
import squidpy as sq
import yaml

# %%
# 步骤：读取配置，设定质控阈值
# 输入：config.yaml
# 输出：cfg、out_dir、阈值
cfg = yaml.safe_load(open("config.yaml"))
out_dir = Path(cfg["out_dir"])
min_counts = 50 if cfg["resolution"] == "cell" else 500    # 单细胞级平台 counts 少，阈值放低
min_cells = 10
max_pct_mt = 20

# %%
# 步骤：读入原始数据，计算质控指标
# 输入：results/01_raw.h5ad
# 输出：adata，obs 里多了 total_counts、n_genes_by_counts、pct_counts_mt
adata = sc.read_h5ad(out_dir / "01_raw.h5ad")
adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, inplace=True)
adata

# %%
# 步骤：按样本看质控分布，再决定阈值
# 输入：adata
# 输出：图 results/qc_violin.png
sc.pl.violin(adata, ["total_counts", "n_genes_by_counts", "pct_counts_mt"], groupby="sample", show=False)
plt.savefig(out_dir / "qc_violin.png", bbox_inches="tight")

# %%
# 步骤：过滤低质量 spot/细胞 和低表达基因
# 输入：adata
# 输出：adata（过滤后），打印过滤后的大小
sc.pp.filter_cells(adata, min_counts=min_counts)
sc.pp.filter_genes(adata, min_cells=min_cells)
adata = adata[adata.obs["pct_counts_mt"] < max_pct_mt].copy()
adata

# %%
# 步骤：标准化和 log 转换
# 输入：adata（X 为 counts）
# 输出：adata（X 为 log1p 标准化值，原始 counts 仍在 layers['counts']）
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# %%
# 步骤：全转录组平台选 2000 个高变基因；靶向 panel（基因数 ≤ 1000）全部基因都用
# 输入：adata
# 输出：adata.var['highly_variable']
if adata.n_vars > 1000:
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")
else:
    adata.var["highly_variable"] = True
adata.var["highly_variable"].sum()

# %%
# 步骤：PCA；多样本时用 Harmony 去批次
# 输入：adata
# 输出：obsm['X_pca']；多样本时还有 obsm['X_pca_harmony']；use_rep（后续用哪个降维结果）
sc.pp.pca(adata, n_comps=30, mask_var="highly_variable", random_state=cfg["seed"])
use_rep = "X_pca"
if adata.obs["sample"].nunique() > 1:
    sc.external.pp.harmony_integrate(adata, "sample")
    use_rep = "X_pca_harmony"
use_rep

# %%
# 步骤：建 KNN 图、UMAP、Leiden 聚类
# 输入：adata、use_rep
# 输出：obs['leiden']、obsm['X_umap']
sc.pp.neighbors(adata, use_rep=use_rep, random_state=cfg["seed"])
sc.tl.umap(adata, random_state=cfg["seed"])
sc.tl.leiden(adata, resolution=cfg["leiden_resolution"], flavor="igraph", n_iterations=2, random_state=cfg["seed"])
adata.obs["leiden"].value_counts()

# %%
# 步骤：画 UMAP（按聚类、按样本）
# 输入：adata
# 输出：图 results/umap_leiden.png
sc.pl.umap(adata, color=["leiden", "sample"], show=False)
plt.savefig(out_dir / "umap_leiden.png", bbox_inches="tight")

# %%
# 步骤：画空间分布（每个样本一张）
# 输入：adata（obsm['spatial_um']）
# 输出：图 results/spatial_leiden.png
sq.pl.spatial_scatter(adata, color="leiden", library_key="sample", spatial_key="spatial_um", shape=None)
plt.savefig(out_dir / "spatial_leiden.png", bbox_inches="tight")

# %%
# 步骤：找每个聚类的 marker 基因
# 输入：adata
# 输出：表 results/markers_leiden.csv
sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
markers = sc.get.rank_genes_groups_df(adata, group=None)
markers.to_csv(out_dir / "markers_leiden.csv", index=False)
markers.groupby("group", observed=True).head(5)

# %%
# 步骤：保存
# 输入：adata
# 输出：results/02_clustered.h5ad
adata.write_h5ad(out_dir / "02_clustered.h5ad")
