# 模块 C · 准备单细胞参考   [Python]
# 对应课程：第4、13、20课
# 作用：按细胞类型下采样，去掉太少的类型，统一列名。RCTD、cell2location、标签转移都从这个文件开始。
# 要求：参考 h5ad 里有原始整数 counts（在 X 或 layers['counts']），并且有细胞类型列。

# %%
# 步骤：导入包
# 输入：无
# 输出：无
from pathlib import Path
import scanpy as sc
import yaml

# %%
# 步骤：读取配置，设定下采样参数
# 输入：config.yaml
# 输出：cfg、out_dir、celltype_key、max_cells / min_cells
cfg = yaml.safe_load(open("config.yaml"))
out_dir = Path(cfg["out_dir"])
celltype_key = cfg["scrna_celltype_key"]
max_cells = 1000     # 每个细胞类型最多保留的细胞数（RCTD 上限是 10000，1000 已经够用）
min_cells = 25       # 少于 25 个细胞的类型 RCTD 会报错，直接去掉

# %%
# 步骤：读入参考，确认 counts 所在位置
# 输入：scrna_reference（h5ad）
# 输出：ref；layers['counts'] 为原始 counts
ref = sc.read_h5ad(cfg["scrna_reference"])
if "counts" not in ref.layers:
    ref.layers["counts"] = ref.X.copy()
ref

# %%
# 步骤：查看每个细胞类型的细胞数
# 输入：ref
# 输出：打印的计数表
ref.obs[celltype_key].value_counts()

# %%
# 步骤：去掉细胞太少的类型，每个类型随机保留最多 max_cells 个细胞
# 输入：ref
# 输出：ref（下采样后）
counts = ref.obs[celltype_key].value_counts()
ref = ref[ref.obs[celltype_key].isin(counts[counts >= min_cells].index)].copy()
keep = ref.obs.sample(frac=1, random_state=cfg["seed"]).groupby(celltype_key, observed=True).head(max_cells).index
ref = ref[keep].copy()
ref.obs[celltype_key].value_counts()

# %%
# 步骤：统一细胞类型列名为 cell_type（RCTD 不接受名字里有 "/"）
# 输入：ref
# 输出：ref.obs['cell_type']
ref.obs["cell_type"] = ref.obs[celltype_key].astype(str).str.replace("/", "_")
ref.var_names_make_unique()

# %%
# 步骤：保存
# 输入：ref
# 输出：results/ref_prepared.h5ad
ref.write_h5ad(out_dir / "ref_prepared.h5ad")
