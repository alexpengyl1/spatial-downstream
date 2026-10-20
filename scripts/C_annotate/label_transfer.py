# 模块 C · 单细胞参考 → 单细胞级空间数据的标签转移   [Python]
# 对应课程：第20课
# 适用：Xenium / CosMx / HD 细胞分割结果 / Stereo-seq cellbin + 有单细胞参考
# 方法：只用两者共有的基因 → 合并 → PCA → Harmony 去平台差异 → KNN（按距离加权）投票预测细胞类型。
#      课程里手写的加权投票循环，这里换成 sklearn 的 KNeighborsClassifier(weights="distance")。

# %%
# 步骤：导入包
# 输入：无
# 输出：无
from pathlib import Path
import matplotlib.pyplot as plt
import scanpy as sc
from sklearn.neighbors import KNeighborsClassifier
import squidpy as sq
import yaml

# %%
# 步骤：读取配置
# 输入：config.yaml
# 输出：cfg、out_dir、n_neighbors
cfg = yaml.safe_load(open("config.yaml"))
out_dir = Path(cfg["out_dir"])
n_neighbors = 30

# %%
# 步骤：读入参考和空间数据，只保留共有基因
# 输入：results/ref_prepared.h5ad、results/02_clustered.h5ad
# 输出：ref、query（基因相同）、shared
ref = sc.read_h5ad(out_dir / "ref_prepared.h5ad")
adata = sc.read_h5ad(out_dir / "02_clustered.h5ad")
shared = adata.var_names.intersection(ref.var_names)
ref = ref[:, shared].copy()
query = adata[:, shared].copy()
len(shared)

# %%
# 步骤：两份数据在共有基因上用同样的方式标准化（从 counts 开始）
# 输入：ref、query
# 输出：ref、query（X 为 log1p 标准化值）
for data in (ref, query):
    data.X = data.layers["counts"].copy()
    sc.pp.normalize_total(data, target_sum=1e4)
    sc.pp.log1p(data)

# %%
# 步骤：合并后做 PCA 和 Harmony（按数据来源去批次）
# 输入：ref、query
# 输出：combined（obsm['X_pca_harmony']）
combined = sc.concat({"ref": ref, "query": query}, label="dataset")
sc.pp.scale(combined, max_value=10)
sc.pp.pca(combined, n_comps=30, random_state=cfg["seed"])
sc.external.pp.harmony_integrate(combined, "dataset")

# %%
# 步骤：用参考细胞训练 KNN 分类器，预测空间细胞的类型和置信度
# 输入：combined
# 输出：adata.obs['cell_type']、adata.obs['cell_type_score']
is_ref = (combined.obs["dataset"] == "ref").to_numpy()
embedding = combined.obsm["X_pca_harmony"]
knn = KNeighborsClassifier(n_neighbors=n_neighbors, weights="distance")
knn.fit(embedding[is_ref], ref.obs["cell_type"])      # combined 里参考细胞在前，顺序和 ref 一致
adata.obs["cell_type"] = knn.predict(embedding[~is_ref])
adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")
adata.obs["cell_type_score"] = knn.predict_proba(embedding[~is_ref]).max(axis=1)
adata.obs["cell_type"].value_counts()

# %%
# 步骤：对照聚类看预测结果是否合理
# 输入：adata
# 输出：打印的聚类 x 细胞类型交叉表
sc.metrics.confusion_matrix("leiden", "cell_type", adata.obs).round(2)

# %%
# 步骤：画细胞类型的空间分布
# 输入：adata
# 输出：图 results/spatial_cell_type.png
sq.pl.spatial_scatter(adata, color="cell_type", library_key="sample", spatial_key="spatial_um", shape=None)
plt.savefig(out_dir / "spatial_cell_type.png", bbox_inches="tight")

# %%
# 步骤：保存
# 输入：adata
# 输出：results/03_annotated.h5ad
adata.write_h5ad(out_dir / "03_annotated.h5ad")
