# 模块 A · 读入 CODEX / IMC / mIHC 单细胞表（空间蛋白组）   [Python]
# 对应课程：第26课
# 输入格式：每个样本一个 CSV，每行一个细胞，包含中心坐标列和每个 marker 的平均强度列（QuPath、HALO、MCMICRO 等都能导出）。

# %%
# 步骤：导入包
# 输入：无
# 输出：无
from pathlib import Path
import pandas as pd
import scanpy as sc
import yaml

# %%
# 步骤：读取配置、样本表，并设定坐标列和 marker 列
# 输入：config.yaml；samples.tsv（path 指向每个样本的 CSV）
# 输出：cfg、samples、out_dir、x_col / y_col / marker_cols
cfg = yaml.safe_load(open("config.yaml"))
out_dir = Path(cfg["out_dir"])
out_dir.mkdir(parents=True, exist_ok=True)
samples = pd.read_csv(cfg["samples"], sep="\t")
x_col, y_col = "Centroid X", "Centroid Y"       # 按自己的 CSV 修改
marker_cols = None                               # None = 除坐标外的所有数值列；也可以写成 ["CD3", "CD20", ...]
samples

# %%
# 步骤：逐个样本读入，构建 AnnData，坐标换算成 µm
# 输入：每个样本的 CSV
# 输出：adatas（样本名 -> AnnData），X 为原始强度
adatas = {}
for row in samples.itertuples():
    table = pd.read_csv(row.path)
    table.index = table.index.astype(str)
    markers = marker_cols or [c for c in table.select_dtypes("number").columns if c not in (x_col, y_col)]
    adata = sc.AnnData(X=table[markers].to_numpy(dtype="float32"), obs=table.drop(columns=markers), var=pd.DataFrame(index=markers))
    adata.obs["condition"] = row.condition
    adata.obsm["spatial"] = table[[x_col, y_col]].to_numpy(dtype="float64")
    adata.obsm["spatial_um"] = adata.obsm["spatial"] * cfg["pixel_size_um"]
    adatas[row.sample] = adata

# %%
# 步骤：合并样本，保存原始强度
# 输入：adatas
# 输出：adata（细胞 x marker）；layers['raw']
adata = sc.concat(adatas, label="sample", index_unique="_")
adata.layers["raw"] = adata.X.copy()
adata

# %%
# 步骤：保存
# 输入：adata
# 输出：results/01_raw.h5ad
adata.write_h5ad(out_dir / "01_raw.h5ad")
