# 模块 G · 空间趋势：以某个区域为中心，基因 / 细胞组成随距离如何变化   [Python，需要 SOAPy_st]
# 对应课程：第10课（课程叫"细胞轨迹"，实际是空间距离上的趋势，不是发育轨迹）
# 适用：所有 RNA 平台；最常见的用法是"肿瘤边界向外，免疫细胞如何变化"
# 注意：SOAPy 的距离单位是 obsm['spatial'] 的单位，这里换成 µm；radius / cut 都按 µm 理解。

# %%
# 步骤：导入包
# 输入：无
# 输出：无
from pathlib import Path
import matplotlib.pyplot as plt
import scanpy as sc
import SOAPy_st as sp
import yaml

# %%
# 步骤：读取配置，指定中心区域
# 输入：config.yaml
# 输出：cfg、out_dir、center_key、center_value、radius_um
cfg = yaml.safe_load(open("config.yaml"))
out_dir = Path(cfg["out_dir"])
center_key = "niche"          # 用哪个分组定义中心区域（niche / domain / cell_type）
center_value = "niche_1"      # 中心区域，比如肿瘤核心所在的 niche
radius_um = 1000              # 从区域边界向外看多远

# %%
# 步骤：读入数据，取一个样本（区域掩膜是按单张切片画的），坐标换成 µm
# 输入：results/05_niche.h5ad
# 输出：adata（单样本）
adata = sc.read_h5ad(out_dir / "05_niche.h5ad")
sample = adata.obs["sample"].unique()[0]
adata = adata[adata.obs["sample"] == sample].copy()
adata.obsm["spatial"] = adata.obsm["spatial_um"]

# %%
# 步骤：根据中心区域生成掩膜
# 输入：adata
# 输出：mask（二维数组）；图 results/tendency_mask.png
mask = sp.tl.get_mask_from_domain(adata, clusters=center_value, KSize=35, cluster_key=center_key)
plt.imshow(mask, cmap="gray")
plt.savefig(out_dir / "tendency_mask.png", bbox_inches="tight")

# %%
# 步骤：哪些基因随离区域的距离显著变化（Spearman 相关）
# 输入：adata、mask
# 输出：表 results/tendency_spearman.csv
spearman = sp.tl.spearman_correlation(adata, mask, radius=radius_um, num=5)
spearman.to_csv(out_dir / "tendency_spearman.csv")
spearman.head(10)

# %%
# 步骤：拟合每个基因随距离的变化曲线，并把曲线相似的基因聚成一类
# 输入：adata、mask
# 输出：图 results/tendency_gene_clusters.png
sp.tl.spatial_tendency(adata, mask, radius=radius_um, location="out", frac=5)
sp.tl.gene_cluster(adata=adata, k=8, range_min=0.03, fdr=True, pvalue=0.05)
sp.pl.show_curves_cluster(adata)
plt.savefig(out_dir / "tendency_gene_clusters.png", bbox_inches="tight")
