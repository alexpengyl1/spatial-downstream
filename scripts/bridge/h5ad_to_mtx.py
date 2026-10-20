# 桥接 · Python → R：把 h5ad 导出成 R 能直接读的文件   [Python → 文件]
# 用在：Python 处理到一半，下一步要用 R 包（BANKSY、Seurat、SpaGene、CellChat、MultiNicheNet 等）
# 导出内容（results/exchange/）：
#   counts.mtx    原始 counts，基因 x 细胞（Matrix Market 格式，R 用 Matrix::readMM 读）
#   genes.txt     基因名，一行一个
#   barcodes.txt  条形码，一行一个
#   obs.csv       细胞元数据（sample、condition、cell_type、niche ……）
#   coords.csv    µm 坐标（x, y）

# %%
# 步骤：导入包
# 输入：无
# 输出：无
from pathlib import Path
import pandas as pd
import scanpy as sc
from scipy import io, sparse
import yaml

# %%
# 步骤：读取配置，指定要导出的文件
# 输入：config.yaml
# 输出：out_dir、input_file、exchange_dir
cfg = yaml.safe_load(open("config.yaml"))
out_dir = Path(cfg["out_dir"])
input_file = out_dir / "03_annotated.h5ad"      # 按需要改成别的阶段文件
exchange_dir = out_dir / "exchange"
exchange_dir.mkdir(parents=True, exist_ok=True)

# %%
# 步骤：读入数据
# 输入：input_file
# 输出：adata
adata = sc.read_h5ad(input_file)
adata

# %%
# 步骤：写出 counts（转置成基因 x 细胞，和 R / Seurat 的习惯一致）、基因名、条形码
# 输入：adata.layers['counts']
# 输出：counts.mtx、genes.txt、barcodes.txt
io.mmwrite(exchange_dir / "counts.mtx", sparse.csr_matrix(adata.layers["counts"]).T)
pd.Series(adata.var_names).to_csv(exchange_dir / "genes.txt", index=False, header=False)
pd.Series(adata.obs_names).to_csv(exchange_dir / "barcodes.txt", index=False, header=False)

# %%
# 步骤：写出元数据和坐标
# 输入：adata.obs、adata.obsm['spatial_um']
# 输出：obs.csv、coords.csv
adata.obs.to_csv(exchange_dir / "obs.csv")
pd.DataFrame(adata.obsm["spatial_um"], index=adata.obs_names, columns=["x", "y"]).to_csv(exchange_dir / "coords.csv")
