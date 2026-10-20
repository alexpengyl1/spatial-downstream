# 模块 F · 只在目标细胞的邻域里做配体-受体分析   [Python]
# 对应课程：第25课
# 适用：单细胞级平台（HD 细胞分割、Xenium、CosMx）
# 思路：先找出目标细胞（比如巨噬细胞）和它半径 radius_um 内的所有邻居，只用这些细胞跑 CellPhoneDB 置换检验。
#      这样得到的是"真正挨在一起的细胞之间"的通讯，而不是整张切片里两类细胞的平均。
# 说明：课程用的是 squidpy 的 sq.gr.ligrec；在 squidpy 1.6.5 上测试时，即使是人为构造的简单例子，p 值也全部是 NaN。
#      所以这里改用 LIANA+ 里实现的 CellPhoneDB 方法，结果的含义相同。

# %%
# 步骤：导入包
# 输入：无
# 输出：无
from pathlib import Path
import liana as li
import numpy as np
import scanpy as sc
import squidpy as sq
import yaml

# %%
# 步骤：读取配置，指定目标细胞类型
# 输入：config.yaml
# 输出：cfg、out_dir、target_cell_type、cluster_key、resource
cfg = yaml.safe_load(open("config.yaml"))
out_dir = Path(cfg["out_dir"])
target_cell_type = "Macrophage"      # 改成自己关心的细胞类型
cluster_key = "cell_type"
resource = "consensus" if cfg["species"] == "human" else "mouseconsensus"

# %%
# 步骤：读入数据，按半径建邻居图
# 输入：results/05_niche.h5ad
# 输出：adata.obsp['spatial_connectivities']
adata = sc.read_h5ad(out_dir / "05_niche.h5ad")
sq.gr.spatial_neighbors(adata, library_key="sample", spatial_key="spatial_um", coord_type="generic", radius=cfg["radius_um"])

# %%
# 步骤：找出目标细胞及其所有邻居，取子集
# 输入：adata
# 输出：nhood（目标细胞 + 邻居），打印其中各类细胞的数量
is_target = (adata.obs[cluster_key] == target_cell_type).to_numpy()
neighbor_index = adata.obsp["spatial_connectivities"][is_target].nonzero()[1]
keep = np.union1d(np.where(is_target)[0], neighbor_index)
nhood = adata[keep].copy()
nhood.obs[cluster_key].value_counts()

# %%
# 步骤：在邻域子集上做 CellPhoneDB 置换检验（配受体库来自 LIANA，首次运行需要联网下载）
# 输入：nhood（X 为 log 标准化表达）
# 输出：nhood.uns['liana_res']（source、target、配体、受体、lr_means、cellphone_pvals）
li.mt.cellphonedb(nhood, groupby=cluster_key, resource_name=resource, expr_prop=0.1, use_raw=False, n_perms=1000, seed=cfg["seed"], verbose=True)
result = nhood.uns["liana_res"]
result.to_csv(out_dir / f"ligrec_{target_cell_type}.csv", index=False)
result.head()

# %%
# 步骤：画目标细胞作为信号发送方的显著配受体对（p ≤ 0.05，按平均表达取前 20 个）
# 输入：nhood
# 输出：图 results/ligrec_<目标细胞>.png
plot = li.pl.dotplot(adata=nhood, colour="lr_means", size="cellphone_pvals", inverse_size=True, source_labels=[target_cell_type], target_labels=list(nhood.obs[cluster_key].cat.categories), filter_fun=lambda x: x["cellphone_pvals"] <= 0.05, top_n=20, orderby="lr_means", orderby_ascending=False, figure_size=(8, 7), return_fig=True)
plot.save(out_dir / f"ligrec_{target_cell_type}.png", verbose=False)
