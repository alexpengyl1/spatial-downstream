# 模块 F · 用 COMMOT 推断空间信号流向（最优传输）   [Python，建议单独的 commot 环境]
# 对应课程：第7课
# 适用：低精度平台（Visium、Stereo-seq bin50）。单细胞级数据点太多，内存会爆，要先下采样或只取一块区域。
# 结果：每个通路的发送 / 接收强度，以及信号流向的箭头图

# %%
# 步骤：导入包
# 输入：无
# 输出：无
from pathlib import Path
import commot as ct
import matplotlib.pyplot as plt
import scanpy as sc
import yaml

# %%
# 步骤：读取配置，设定参数
# 输入：config.yaml
# 输出：cfg、out_dir、distance_um、n_pathways
cfg = yaml.safe_load(open("config.yaml"))
out_dir = Path(cfg["out_dir"])
distance_um = 200      # 信号最远能传多远（µm）；分泌型信号一般 100–250µm
n_pathways = 6         # 画信号最强的前几个通路

# %%
# 步骤：读入数据，取一个样本（COMMOT 按坐标距离计算，必须单样本运行），坐标换成 µm
# 输入：results/05_niche.h5ad
# 输出：adata（单样本，obsm['spatial'] 为 µm）
adata = sc.read_h5ad(out_dir / "05_niche.h5ad")
sample = adata.obs["sample"].unique()[0]          # 换样本就改这里
adata = adata[adata.obs["sample"] == sample].copy()
adata.obsm["spatial"] = adata.obsm["spatial_um"]
adata

# %%
# 步骤：取 CellChat 的分泌型配受体库，只保留在 ≥5% spot 中表达的配受体
# 输入：adata
# 输出：lr_table（配体、受体、通路）
lr_table = ct.pp.ligand_receptor_database(species=cfg["species"], signaling_type="Secreted Signaling", database="CellChat")
lr_table = ct.pp.filter_lr_database(lr_table, adata, min_cell_pct=0.05)
lr_table.shape

# %%
# 步骤：计算空间通讯（这一步最慢）
# 输入：adata、lr_table
# 输出：adata.obsm['commot-cellchat-sum-sender'] / ['commot-cellchat-sum-receiver'] 等
ct.tl.spatial_communication(adata, database_name="cellchat", df_ligrec=lr_table, dis_thr=distance_um, heteromeric=True, pathway_sum=True)

# %%
# 步骤：按总发送量给通路排序
# 输入：adata.obsm['commot-cellchat-sum-sender']
# 输出：top_pathways（列表）
sender = adata.obsm["commot-cellchat-sum-sender"]          # 列名为 s-<通路> 或 s-<配体-受体>
pathway_cols = [f"s-{p}" for p in lr_table.iloc[:, 2].unique() if f"s-{p}" in sender.columns]
top_pathways = sender[pathway_cols].sum().nlargest(n_pathways).index.str[2:].tolist()
top_pathways

# %%
# 步骤：对每个通路计算信号方向，画流向图（背景按 niche 上色）
# 输入：adata、top_pathways
# 输出：每个通路一张图 results/commot_<通路>.png
for pathway in top_pathways:
    ct.tl.communication_direction(adata, database_name="cellchat", pathway_name=pathway, k=5)
    ct.pl.plot_cell_communication(adata, database_name="cellchat", pathway_name=pathway, plot_method="grid", background="cluster", clustering="niche", normalize_v=True, normalize_v_quantile=0.995, grid_density=0.4, scale=0.00001, ndsize=8, filename=str(out_dir / f"commot_{pathway}.png"))
    plt.close("all")
