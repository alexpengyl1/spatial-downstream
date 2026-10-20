# 模块 A · 读入 NanoString CosMx（单细胞级，靶向 panel）   [Python]
# 对应课程：第19课
# 说明：用 CenterX/Y_global_px（全片坐标），不要用 local_px：不同 FOV 的 local 坐标会重叠。

# %%
# 步骤：导入包
# 输入：无
# 输出：无
from pathlib import Path
import pandas as pd
import scanpy as sc
from scipy import sparse
import yaml

# %%
# 步骤：读取配置和样本表
# 输入：config.yaml；samples.tsv（path 指向含 *exprMat_file.csv 和 *metadata_file.csv 的目录）
# 输出：cfg、samples、out_dir
cfg = yaml.safe_load(open("config.yaml"))
out_dir = Path(cfg["out_dir"])
out_dir.mkdir(parents=True, exist_ok=True)
samples = pd.read_csv(cfg["samples"], sep="\t")
samples

# %%
# 步骤：逐个样本读入表达矩阵和元数据，去掉背景（cell_ID = 0），用"FOV_细胞号"对齐两张表
# 输入：*exprMat_file.csv、*metadata_file.csv
# 输出：tables（样本名 -> (表达表, 元数据表)）
tables = {}
for row in samples.itertuples():
    expr = pd.read_csv(next(Path(row.path).glob("*exprMat_file.csv")))
    meta = pd.read_csv(next(Path(row.path).glob("*metadata_file.csv")))
    expr = expr[expr["cell_ID"] != 0]
    expr.index = expr["fov"].astype(str) + "_" + expr["cell_ID"].astype(str)
    meta.index = meta["fov"].astype(str) + "_" + meta["cell_ID"].astype(str)
    tables[row.sample] = (expr, meta.loc[expr.index])

# %%
# 步骤：每个样本构建 AnnData（去掉阴性探针列），坐标换算成 µm
# 输入：tables
# 输出：adatas（样本名 -> AnnData），obsm['spatial']（像素）、obsm['spatial_um']（µm）
adatas = {}
for row in samples.itertuples():
    expr, meta = tables[row.sample]
    genes = [g for g in expr.columns[2:] if not g.startswith(("NegPrb", "Negative", "SystemControl"))]
    adata = sc.AnnData(X=sparse.csr_matrix(expr[genes].to_numpy()), obs=meta, var=pd.DataFrame(index=genes))
    adata.obs["condition"] = row.condition
    adata.obsm["spatial"] = adata.obs[["CenterX_global_px", "CenterY_global_px"]].to_numpy()
    adata.obsm["spatial_um"] = adata.obsm["spatial"] * cfg["pixel_size_um"]
    adatas[row.sample] = adata

# %%
# 步骤：合并样本，保存原始 counts
# 输入：adatas
# 输出：adata；layers['counts']
adata = sc.concat(adatas, label="sample", index_unique="_")
adata.layers["counts"] = adata.X.copy()
adata

# %%
# 步骤：保存
# 输入：adata
# 输出：results/01_raw.h5ad
adata.write_h5ad(out_dir / "01_raw.h5ad")
