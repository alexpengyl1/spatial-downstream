# 模块 C · 用 cell2location 做去卷积（RCTD 的 Python 替代方案）   [Python，需要 cell2location 环境，最好有 GPU]
# 对应课程：第4、5课
# 适用：低精度平台（Visium、Stereo-seq bin50）+ 有单细胞参考
# 结果：obsm['cell_props']（细胞类型比例）、obsm['cell_abundance']（每个 spot 每类细胞的估计个数）、obs['cell_type']
# 注意：spot 训练 30000 轮在 CPU 上要几个小时；先用 max_epochs 小的值把流程跑通。

# %%
# 步骤：导入包
# 输入：无
# 输出：无
from pathlib import Path
from cell2location.models import Cell2location, RegressionModel
from cell2location.utils.filtering import filter_genes
import matplotlib.pyplot as plt
import scanpy as sc
import squidpy as sq
import yaml

# %%
# 步骤：读取配置，设定模型参数
# 输入：config.yaml
# 输出：cfg、out_dir、参数
cfg = yaml.safe_load(open("config.yaml"))
out_dir = Path(cfg["out_dir"])
n_cells_per_spot = 10        # 每个 spot 预期的细胞数：Visium 一般 5–20，按组织密度调整
ref_epochs = 250
spatial_epochs = 30000

# %%
# 步骤：读入参考和空间数据，两者都把 X 换回原始 counts，空间数据去掉线粒体基因
# 输入：results/ref_prepared.h5ad、results/02_clustered.h5ad
# 输出：ref、spatial
ref = sc.read_h5ad(out_dir / "ref_prepared.h5ad")
ref.X = ref.layers["counts"].copy()
spatial = sc.read_h5ad(out_dir / "02_clustered.h5ad")
spatial.X = spatial.layers["counts"].copy()
spatial = spatial[:, ~spatial.var_names.str.upper().str.startswith("MT-")].copy()

# %%
# 步骤：参考数据去掉信息量少的基因
# 输入：ref
# 输出：ref（基因更少）
selected = filter_genes(ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
ref = ref[:, selected].copy()

# %%
# 步骤：用负二项回归估计每个细胞类型的表达特征
# 输入：ref
# 输出：signatures（基因 x 细胞类型）
RegressionModel.setup_anndata(ref, batch_key=cfg["scrna_batch_key"], labels_key="cell_type")
ref_model = RegressionModel(ref)
ref_model.train(max_epochs=ref_epochs)
ref = ref_model.export_posterior(ref, sample_kwargs={"num_samples": 1000, "batch_size": 2500})
signatures = ref.varm["means_per_cluster_mu_fg"].copy()
signatures.columns = ref.uns["mod"]["factor_names"]

# %%
# 步骤：取空间数据和表达特征共有的基因
# 输入：spatial、signatures
# 输出：spatial、signatures（同样的基因）
shared = spatial.var_names.intersection(signatures.index)
spatial = spatial[:, shared].copy()
signatures = signatures.loc[shared]
len(shared)

# %%
# 步骤：训练 cell2location 空间模型（这一步最慢）
# 输入：spatial、signatures
# 输出：spatial_model（训练好的模型）
Cell2location.setup_anndata(spatial, batch_key="sample")
spatial_model = Cell2location(spatial, cell_state_df=signatures, N_cells_per_location=n_cells_per_spot, detection_alpha=20)
spatial_model.train(max_epochs=spatial_epochs, batch_size=None, train_size=1)

# %%
# 步骤：导出每个 spot 的细胞丰度，换算成比例，放回聚类后的数据
# 输入：spatial_model
# 输出：adata（obsm['cell_abundance']、obsm['cell_props']、obs['cell_type']）
spatial = spatial_model.export_posterior(spatial, sample_kwargs={"num_samples": 1000, "batch_size": spatial_model.adata.n_obs})
abundance = spatial.obsm["q05_cell_abundance_w_sf"]
abundance.columns = spatial.uns["mod"]["factor_names"]
adata = sc.read_h5ad(out_dir / "02_clustered.h5ad")
adata.obsm["cell_abundance"] = abundance.loc[adata.obs_names]
adata.obsm["cell_props"] = adata.obsm["cell_abundance"].div(adata.obsm["cell_abundance"].sum(axis=1), axis=0)
adata.obs["cell_type"] = adata.obsm["cell_props"].idxmax(axis=1).astype("category")

# %%
# 步骤：画每类细胞丰度的空间分布
# 输入：adata
# 输出：图 results/c2l_abundance.png
for cell_type in adata.obsm["cell_abundance"].columns:
    adata.obs[f"abund_{cell_type}"] = adata.obsm["cell_abundance"][cell_type]
sq.pl.spatial_scatter(adata, color=[f"abund_{c}" for c in adata.obsm["cell_abundance"].columns], library_key="sample", spatial_key="spatial_um", shape=None)
plt.savefig(out_dir / "c2l_abundance.png", bbox_inches="tight")

# %%
# 步骤：保存
# 输入：adata
# 输出：results/03_annotated.h5ad
adata.write_h5ad(out_dir / "03_annotated.h5ad")
