# 模块 C · 第 1 步（共 3 步）：把空间数据和单细胞参考导出给 R 的 RCTD   [Python → 文件]
# 对应课程：第4、5、13课
# 联合流程：export_for_rctd.py（Python 导出）→ run_rctd.Rmd（R 运行 RCTD）→ import_rctd.py（Python 读回结果）
# 交换格式：counts 用 Matrix Market（.mtx，基因 x 细胞），基因名和条形码各存一个 txt，坐标和细胞类型存 csv。
#           这种格式 R 和 Python 都能直接读，不依赖 SeuratDisk 或 zellkonverter。

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
# 步骤：读取配置，建立交换目录
# 输入：config.yaml
# 输出：cfg、out_dir、rctd_dir（results/rctd）
cfg = yaml.safe_load(open("config.yaml"))
out_dir = Path(cfg["out_dir"])
rctd_dir = out_dir / "rctd"
(rctd_dir / "reference").mkdir(parents=True, exist_ok=True)

# %%
# 步骤：导出单细胞参考（counts、基因、条形码、细胞类型）
# 输入：results/ref_prepared.h5ad
# 输出：results/rctd/reference/ 下的 counts.mtx、genes.txt、barcodes.txt、cell_types.csv
ref = sc.read_h5ad(out_dir / "ref_prepared.h5ad")
io.mmwrite(rctd_dir / "reference" / "counts.mtx", sparse.csr_matrix(ref.layers["counts"]).T)
pd.Series(ref.var_names).to_csv(rctd_dir / "reference" / "genes.txt", index=False, header=False)
pd.Series(ref.obs_names).to_csv(rctd_dir / "reference" / "barcodes.txt", index=False, header=False)
ref.obs[["cell_type"]].to_csv(rctd_dir / "reference" / "cell_types.csv")

# %%
# 步骤：按样本导出空间数据（counts、基因、条形码、µm 坐标）；RCTD 每个样本单独运行
# 输入：results/02_clustered.h5ad
# 输出：results/rctd/spatial/<样本>/ 下的 counts.mtx、genes.txt、barcodes.txt、coords.csv
adata = sc.read_h5ad(out_dir / "02_clustered.h5ad")
for sample in adata.obs["sample"].unique():
    sub = adata[adata.obs["sample"] == sample]
    sample_dir = rctd_dir / "spatial" / sample
    sample_dir.mkdir(parents=True, exist_ok=True)
    io.mmwrite(sample_dir / "counts.mtx", sparse.csr_matrix(sub.layers["counts"]).T)
    pd.Series(sub.var_names).to_csv(sample_dir / "genes.txt", index=False, header=False)
    pd.Series(sub.obs_names).to_csv(sample_dir / "barcodes.txt", index=False, header=False)
    pd.DataFrame(sub.obsm["spatial_um"], index=sub.obs_names, columns=["x", "y"]).to_csv(sample_dir / "coords.csv")
