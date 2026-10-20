# 桥接 · R → Python：把 R 算出的逐细胞结果（CSV，第一列为条形码）合并回 h5ad   [文件 → Python]
# 用在：BANKSY 的 domain、Seurat 的注释、SpaGene 的 pattern 等，任何"一行一个细胞"的结果表
# 原则：按条形码对齐（.reindex），不按行号对齐；R 端如果丢了部分细胞，对应位置是 NaN。

# %%
# 步骤：导入包
# 输入：无
# 输出：无
from pathlib import Path
import pandas as pd
import scanpy as sc
import yaml

# %%
# 步骤：读取配置，指定输入、结果表和输出
# 输入：config.yaml
# 输出：input_file、result_csv、output_file
cfg = yaml.safe_load(open("config.yaml"))
out_dir = Path(cfg["out_dir"])
input_file = out_dir / "03_annotated.h5ad"
result_csv = out_dir / "exchange" / "banksy_domain.csv"
output_file = out_dir / "04_domain.h5ad"

# %%
# 步骤：读入 h5ad 和 R 的结果表
# 输入：input_file、result_csv
# 输出：adata、result（条形码为索引）
adata = sc.read_h5ad(input_file)
result = pd.read_csv(result_csv, index_col=0)
result.head()

# %%
# 步骤：按条形码把结果列合并到 obs；整数和文字列（聚类编号、细胞类型）转成分类变量，小数列保持数值
# 输入：adata、result
# 输出：adata.obs 里多了 result 的各列
for column in result.columns:
    values = result[column]
    if not pd.api.types.is_float_dtype(values):
        values = values.astype(str).astype("category")
    adata.obs[column] = values.reindex(adata.obs_names)
adata.obs[result.columns].head()

# %%
# 步骤：保存
# 输入：adata
# 输出：output_file
adata.write_h5ad(output_file)
