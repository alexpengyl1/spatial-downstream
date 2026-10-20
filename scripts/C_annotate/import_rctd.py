# 模块 C · 第 3 步（共 3 步）：把 R 的 RCTD 结果读回 AnnData   [文件 → Python]
# 对应课程：第4、5、13课
# 结果：obsm['cell_props']（细胞类型比例，低精度平台后续的 niche / 共定位 / 通讯都用它），
#       obs['cell_type']（主导细胞类型；doublet 模式用 RCTD 的 first_type）
# 注意：RCTD 会丢掉 UMI 太少的 spot，这里只保留有结果的 spot。

# %%
# 步骤：导入包
# 输入：无
# 输出：无
from pathlib import Path
import pandas as pd
import scanpy as sc
import yaml

# %%
# 步骤：读取配置和聚类后的数据
# 输入：config.yaml、results/02_clustered.h5ad
# 输出：cfg、out_dir、adata
cfg = yaml.safe_load(open("config.yaml"))
out_dir = Path(cfg["out_dir"])
adata = sc.read_h5ad(out_dir / "02_clustered.h5ad")

# %%
# 步骤：读入所有样本的细胞类型比例，按条形码对齐
# 输入：results/rctd/spatial/*/weights.csv
# 输出：weights（spot x 细胞类型），adata（只保留有 RCTD 结果的 spot）
weights = pd.concat([pd.read_csv(f, index_col=0) for f in (out_dir / "rctd" / "spatial").glob("*/weights.csv")])
adata = adata[adata.obs_names.isin(weights.index)].copy()
adata.obsm["cell_props"] = weights.loc[adata.obs_names]
adata

# %%
# 步骤：给每个 spot 一个主导细胞类型标签
# 输入：weights；doublet 模式下的 results_df.csv
# 输出：obs['cell_type']
adata.obs["cell_type"] = adata.obsm["cell_props"].idxmax(axis=1).astype("category")
if cfg["rctd_mode"] == "doublet":
    results = pd.concat([pd.read_csv(f, index_col=0) for f in (out_dir / "rctd" / "spatial").glob("*/results_df.csv")])
    adata.obs["cell_type"] = results.loc[adata.obs_names, "first_type"].astype("category")
adata.obs["cell_type"].value_counts()

# %%
# 步骤：保存
# 输入：adata
# 输出：results/03_annotated.h5ad
adata.write_h5ad(out_dir / "03_annotated.h5ad")
