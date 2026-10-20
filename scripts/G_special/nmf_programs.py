# 模块 G · NMF 表达程序（把表达拆成若干个"基因程序"，看每个程序的空间分布）   [Python]
# 对应课程：第18课
# 适用：全转录组平台（Visium、HD、Stereo-seq）；靶向 panel 也能跑，但程序数要少。
# 注意：NMF 要求非负输入，用 log 标准化后的值（不能用 scale 之后的值）；HD 数据点多时先下采样拟合，再投影全部数据。

# %%
# 步骤：导入包
# 输入：无
# 输出：无
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
from sklearn.decomposition import NMF
import squidpy as sq
import yaml

# %%
# 步骤：读取配置，设定程序个数
# 输入：config.yaml
# 输出：cfg、out_dir、n_programs、n_fit
cfg = yaml.safe_load(open("config.yaml"))
out_dir = Path(cfg["out_dir"])
n_programs = 10      # 可以试 6 / 10 / 15，看程序是否可解释
n_fit = 50000        # 拟合时最多用多少个 spot / 细胞

# %%
# 步骤：读入数据，只用高变基因
# 输入：results/05_niche.h5ad
# 输出：expr（spot x 高变基因，log 标准化值）
adata = sc.read_h5ad(out_dir / "05_niche.h5ad")
expr = adata[:, adata.var["highly_variable"]].copy()
expr.shape

# %%
# 步骤：在下采样的数据上拟合 NMF，再把全部数据投影上去
# 输入：expr
# 输出：usage（spot x 程序）、loadings（程序 x 基因）
fit_cells = sc.pp.subsample(expr, n_obs=min(n_fit, expr.n_obs), random_state=cfg["seed"], copy=True)
nmf = NMF(n_components=n_programs, init="nndsvda", max_iter=500, random_state=cfg["seed"])
nmf.fit(fit_cells.X)
programs = [f"NMF_{i + 1}" for i in range(n_programs)]
usage = pd.DataFrame(nmf.transform(expr.X), index=adata.obs_names, columns=programs)
loadings = pd.DataFrame(nmf.components_, index=programs, columns=expr.var_names)

# %%
# 步骤：保存每个程序权重最高的 20 个基因
# 输入：loadings
# 输出：表 results/nmf_top_genes.csv
top_genes = pd.DataFrame({p: loadings.loc[p].nlargest(20).index for p in programs})
top_genes.to_csv(out_dir / "nmf_top_genes.csv", index=False)
top_genes.head(10)

# %%
# 步骤：画每个程序的空间分布（按行归一化后的使用比例）
# 输入：usage
# 输出：图 results/spatial_nmf.png
adata.obs[programs] = usage.div(usage.sum(axis=1), axis=0).fillna(0)
sq.pl.spatial_scatter(adata, color=programs, library_key="sample", spatial_key="spatial_um", shape=None, cmap="magma")
plt.savefig(out_dir / "spatial_nmf.png", bbox_inches="tight")

# %%
# 步骤：每个 niche 的程序使用量，帮助解释 niche 的功能
# 输入：adata
# 输出：图 results/nmf_by_niche.png
sc.pl.matrixplot(adata, programs, groupby="niche", standard_scale="var", cmap="Reds", show=False)
plt.savefig(out_dir / "nmf_by_niche.png", bbox_inches="tight")
