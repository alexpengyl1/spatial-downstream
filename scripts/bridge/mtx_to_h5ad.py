# 桥接 · 文件 → Python：把 seurat_to_mtx.Rmd（或任何 R 脚本）导出的文件组装成 h5ad   [文件 → Python]
# 用在：Seurat 单细胞参考转成 h5ad 后，交给 prepare_reference.py

# %%
# 步骤：导入包
# 输入：无
# 输出：无
from pathlib import Path
import pandas as pd
import scanpy as sc
from scipy import io, sparse

# %%
# 步骤：指定输入目录和输出文件
# 输入：无
# 输出：exchange_dir、output_file
exchange_dir = Path("results/exchange_seurat")
output_file = Path("ref/scrna.h5ad")

# %%
# 步骤：读入 counts（基因 x 细胞）并转置成 细胞 x 基因，加上基因名、条形码、元数据
# 输入：counts.mtx、genes.txt、barcodes.txt、obs.csv
# 输出：adata（layers['counts'] 为原始 counts）
counts = sparse.csr_matrix(io.mmread(exchange_dir / "counts.mtx").T)
genes = pd.read_csv(exchange_dir / "genes.txt", header=None)[0].to_numpy()
barcodes = pd.read_csv(exchange_dir / "barcodes.txt", header=None)[0].to_numpy()
obs = pd.read_csv(exchange_dir / "obs.csv", index_col=0).loc[barcodes]
adata = sc.AnnData(X=counts, obs=obs, var=pd.DataFrame(index=genes))
adata.layers["counts"] = adata.X.copy()
adata

# %%
# 步骤：保存
# 输入：adata
# 输出：output_file
output_file.parent.mkdir(parents=True, exist_ok=True)
adata.write_h5ad(output_file)
