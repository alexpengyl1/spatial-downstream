# 模块 F · 空间约束的配体-受体共表达（局部 + 全局）   [Python]
# 对应课程：第7、17课（课程用 SpaGene_LR / squidpy ligrec；这里用 LIANA+ 的 bivariate，真正利用了空间距离）
# 适用：所有 RNA 平台。靶向 panel（Xenium / CosMx）先看输出里剩下多少配受体对，太少就没有意义。
# 注意：squidpy 的 sq.gr.ligrec 不使用空间坐标（只按细胞类型做置换检验），课程第17课先建了半径图，但 ligrec 并没有用到它。
#      需要按空间来约束时，用本脚本，或用 ligrec_neighborhood.py 先圈出邻域。

# %%
# 步骤：导入包
# 输入：无
# 输出：无
from pathlib import Path
import liana as li
import matplotlib.pyplot as plt
import scanpy as sc
import squidpy as sq
import yaml

# %%
# 步骤：读取配置，设定空间核参数
# 输入：config.yaml
# 输出：cfg、out_dir、bandwidth_um、resource
cfg = yaml.safe_load(open("config.yaml"))
out_dir = Path(cfg["out_dir"])
bandwidth_um = 200 if cfg["resolution"] == "spot" else 50     # 信号作用范围：低精度 ~ 2 个 spot，单细胞级 ~ 几个细胞
resource = "consensus" if cfg["species"] == "human" else "mouseconsensus"

# %%
# 步骤：读入数据
# 输入：results/05_niche.h5ad
# 输出：adata（X 为 log 标准化表达）
adata = sc.read_h5ad(out_dir / "05_niche.h5ad")
adata.obsm["spatial"] = adata.obsm["spatial_um"]

# %%
# 步骤：每个样本单独计算（LIANA 的邻居图不区分样本），保存全局 Moran's I 排名和局部得分
# 输入：adata
# 输出：results/lr_global_<样本>.csv；results/lr_local_<样本>.h5ad（spot x 配受体对）
for sample in adata.obs["sample"].unique():
    sub = adata[adata.obs["sample"] == sample].copy()
    li.ut.spatial_neighbors(sub, bandwidth=bandwidth_um, cutoff=0.1, kernel="gaussian", set_diag=True)
    lrdata = li.mt.bivariate(sub, resource_name=resource, local_name="cosine", global_name="morans", n_perms=100, nz_prop=0.2, use_raw=False)
    lrdata.var.sort_values("morans", ascending=False).to_csv(out_dir / f"lr_global_{sample}.csv")
    lrdata.write_h5ad(out_dir / f"lr_local_{sample}.h5ad")

# %%
# 步骤：画第一个样本中空间相关最强的 4 个配受体对的局部得分
# 输入：results/lr_local_<第一个样本>.h5ad
# 输出：图 results/lr_local_top.png
sample = adata.obs["sample"].unique()[0]
lrdata = sc.read_h5ad(out_dir / f"lr_local_{sample}.h5ad")
top_pairs = lrdata.var.sort_values("morans", ascending=False).index[:4].tolist()
sq.pl.spatial_scatter(lrdata, color=top_pairs, shape=None, cmap="coolwarm")
plt.savefig(out_dir / "lr_local_top.png", bbox_inches="tight")
