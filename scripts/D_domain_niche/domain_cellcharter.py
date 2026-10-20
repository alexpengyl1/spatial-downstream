# 模块 D · 空间域（结合邻居表达的聚类）   [Python]
# 对应课程：第12课（课程用 BANKSY；R 版 BANKSY 见 domain_banksy.Rmd）
# 适用：所有 RNA 平台；空间蛋白组也能用（use_rep 改成 "X"）
# 思路：把每个细胞自己的降维结果，与 1~n 层空间邻居的平均降维结果拼在一起，再聚类。
#      这和 BANKSY 用 λ 混合"自身 + 邻居"表达是同一个思路。

# %%
# 步骤：导入包
# 输入：无
# 输出：无
from pathlib import Path
import cellcharter as cc
import matplotlib.pyplot as plt
import scanpy as sc
import squidpy as sq
import yaml

# %%
# 步骤：读取配置，设定参数
# 输入：config.yaml
# 输出：cfg、out_dir、n_layers、domain_resolution
cfg = yaml.safe_load(open("config.yaml"))
out_dir = Path(cfg["out_dir"])
n_layers = 3               # 聚合几层邻居：层数越多，空间域越大越平滑
domain_resolution = 0.5

# %%
# 步骤：读入注释后的数据，按样本建 Delaunay 空间邻居图（样本之间不连边），去掉过长的连边
# 输入：results/03_annotated.h5ad
# 输出：adata.obsp['spatial_connectivities']
adata = sc.read_h5ad(out_dir / "03_annotated.h5ad")
sq.gr.spatial_neighbors(adata, library_key="sample", spatial_key="spatial_um", coord_type="generic", delaunay=True)
cc.gr.remove_long_links(adata)

# %%
# 步骤：聚合邻居的降维结果
# 输入：adata（多样本用 X_pca_harmony，单样本用 X_pca；蛋白组没有 PCA，用 None 表示直接用 X）
# 输出：obsm['X_cellcharter']
use_rep = None
if "X_pca_harmony" in adata.obsm:
    use_rep = "X_pca_harmony"
elif "X_pca" in adata.obsm:
    use_rep = "X_pca"
cc.gr.aggregate_neighbors(adata, n_layers=n_layers, use_rep=use_rep, out_key="X_cellcharter")

# %%
# 步骤：在聚合后的表示上聚类，得到空间域
# 输入：adata.obsm['X_cellcharter']
# 输出：obs['domain']
sc.pp.neighbors(adata, use_rep="X_cellcharter", key_added="domain_nn", random_state=cfg["seed"])
sc.tl.leiden(adata, neighbors_key="domain_nn", key_added="domain", resolution=domain_resolution, flavor="igraph", n_iterations=2, random_state=cfg["seed"])
adata.obs["domain"].value_counts()

# %%
# 步骤：画空间域
# 输入：adata
# 输出：图 results/spatial_domain.png
sq.pl.spatial_scatter(adata, color="domain", library_key="sample", spatial_key="spatial_um", shape=None)
plt.savefig(out_dir / "spatial_domain.png", bbox_inches="tight")

# %%
# 步骤：保存
# 输入：adata
# 输出：results/04_domain.h5ad
adata.write_h5ad(out_dir / "04_domain.h5ad")
