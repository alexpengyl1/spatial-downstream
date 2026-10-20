# 模块 C · 没有单细胞参考时，用 marker 基因（或抗体）注释   [Python]
# 对应课程：第2、19、26课
# 单细胞级平台（resolution: cell）：先给每个细胞打分，再按 Leiden 聚类取平均分最高的类型，整簇注释成这个类型。
# 低精度平台（resolution: spot）：spot 是混合的，不硬性注释。打分矩阵存进 obsm['cell_scores']，
#   后续 niche / MISTy 把它当成细胞组成的近似（config 里设 composition_key: cell_scores）；
#   obs['cell_type'] 只记录得分最高的类型，供画图参考。

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
# 步骤：读取配置、marker 表和聚类后的数据
# 输入：config.yaml、markers.yaml、results/02_clustered.h5ad
# 输出：cfg、markers（细胞类型 -> 基因列表）、adata
cfg = yaml.safe_load(open("config.yaml"))
out_dir = Path(cfg["out_dir"])
markers = yaml.safe_load(open(cfg["markers"]))
adata = sc.read_h5ad(out_dir / "02_clustered.h5ad")

# %%
# 步骤：只保留数据里存在的 marker，看每个类型还剩几个
# 输入：markers、adata
# 输出：markers（过滤后），打印每个类型的 marker 数
markers = {cell_type: [g for g in genes if g in adata.var_names] for cell_type, genes in markers.items()}
{cell_type: len(genes) for cell_type, genes in markers.items()}

# %%
# 步骤：每个细胞类型打一个 marker 分
#       靶向 panel 和蛋白 marker 少：对照基因数和分箱数调小，并允许 marker 本身进入对照池（ctrl_as_ref=False），否则会找不到对照基因
# 输入：adata、markers
# 输出：obsm['cell_scores']（spot/细胞 x 细胞类型）
small_panel = adata.n_vars <= 1000
ctrl_size = min(50, adata.n_vars // 4)
n_bins = 5 if small_panel else 25
for cell_type, genes in markers.items():
    sc.tl.score_genes(adata, genes, score_name=f"score_{cell_type}", ctrl_size=ctrl_size, n_bins=n_bins, ctrl_as_ref=not small_panel, random_state=cfg["seed"])
scores = adata.obs[[f"score_{c}" for c in markers]]
scores.columns = list(markers)
adata.obsm["cell_scores"] = scores

# %%
# 步骤：单细胞级平台按聚类注释；低精度平台按 spot 取得分最高的类型
# 输入：adata
# 输出：obs['cell_type']
if cfg["resolution"] == "cell":
    cluster_label = scores.groupby(adata.obs["leiden"], observed=True).mean().idxmax(axis=1)
    adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_label)
else:
    adata.obs["cell_type"] = scores.idxmax(axis=1)
adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")
pd.crosstab(adata.obs["leiden"], adata.obs["cell_type"])

# %%
# 步骤：画 marker 分数，检查注释
# 输入：adata
# 输出：图 results/marker_scores.png
sc.pl.dotplot(adata, markers, groupby="leiden", show=False)
plt.savefig(out_dir / "marker_scores.png", bbox_inches="tight")

# %%
# 步骤：画细胞类型的空间分布
# 输入：adata
# 输出：图 results/spatial_cell_type.png
sq.pl.spatial_scatter(adata, color="cell_type", library_key="sample", spatial_key="spatial_um", shape=None)
plt.savefig(out_dir / "spatial_cell_type.png", bbox_inches="tight")

# %%
# 步骤：保存
# 输入：adata
# 输出：results/03_annotated.h5ad
adata.write_h5ad(out_dir / "03_annotated.h5ad")
