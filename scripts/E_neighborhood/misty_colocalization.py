# 模块 E · 用 MISTy 做细胞类型共定位   [Python]
# 对应课程：第6课（共定位分析，单样本 / 多样本比较）
# 适用：低精度平台的细胞组成（RCTD / cell2location 比例，或 marker 打分）；也可以用在单细胞级平台的邻域组成上。
# 思路：用 spot 自己（intra）、紧邻（juxta）和周围一定范围（para）的细胞组成，
#       通过随机森林预测每一类细胞的丰度；某类细胞对预测贡献越大，说明两者在空间上越相关。
# 多样本：逐样本建模，再汇总。MISTy 只看坐标，每个样本单独建模才能保证邻居都在同一张切片；
#         每个样本的重要性都保存下来，既可以取平均看共同规律，也可以按分组比较（课程第6课的"多样本比较"，
#         以及课件引用的心肌梗死空间多组学文章都是这样做的）。

# %%
# 步骤：导入包
# 输入：无
# 输出：无
from pathlib import Path
import liana as li
from liana.method import genericMistyData
from liana.method.sp import RandomForestModel
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns
import yaml

# %%
# 步骤：读取配置，设定 MISTy 参数，读入数据
# 输入：config.yaml、results/05_niche.h5ad
# 输出：adata、key（用哪个细胞组成）、bandwidth_um、n_neighs
cfg = yaml.safe_load(open("config.yaml"))
out_dir = Path(cfg["out_dir"])
bandwidth_um = 200      # para 视图的范围（µm）：Visium 相邻 spot 间距约 100µm
n_neighs = 6            # juxta 视图用几个最近邻：Visium 六边形排列，正好 6 个
adata = sc.read_h5ad(out_dir / "05_niche.h5ad")
key = cfg["composition_key"] if cfg["resolution"] == "spot" else "nhood_composition"

# %%
# 步骤：每个样本单独构建 MISTy 的三个视图并训练随机森林
# 输入：adata（obsm[key]、obsm['spatial_um']）
# 输出：interactions、metrics（所有样本的结果，带 sample 列）
interactions, metrics = [], []
for sample in adata.obs["sample"].unique():
    comps = li.ut.obsm_to_adata(adata[adata.obs["sample"] == sample], key)
    comps.obsm["spatial"] = adata[adata.obs["sample"] == sample].obsm["spatial_um"]
    misty = genericMistyData(intra=comps, extra=comps, cutoff=0.05, bandwidth=bandwidth_um, n_neighs=n_neighs)
    misty(model=RandomForestModel, n_jobs=cfg["threads"], seed=cfg["seed"])
    interactions.append(misty.uns["interactions"].assign(sample=sample))
    metrics.append(misty.uns["target_metrics"].assign(sample=sample))

# %%
# 步骤：保存每个样本的结果（用于按分组比较）
# 输入：interactions、metrics
# 输出：表 results/misty_interactions_by_sample.csv、results/misty_target_metrics_by_sample.csv
interactions = pd.concat(interactions, ignore_index=True)
metrics = pd.concat(metrics, ignore_index=True)
interactions.to_csv(out_dir / "misty_interactions_by_sample.csv", index=False)
metrics.to_csv(out_dir / "misty_target_metrics_by_sample.csv", index=False)
metrics.groupby("target")[["intra_R2", "multi_R2"]].mean().round(2)

# %%
# 步骤：紧邻视图（juxta）的重要性在各样本间取平均，画成 目标细胞 x 预测细胞 的热图
# 输入：interactions
# 输出：表 results/misty_juxta_mean.csv；图 results/misty_juxta_mean.png
juxta = interactions[interactions["view"] == "juxta"]
juxta_mean = juxta.pivot_table(index="target", columns="predictor", values="importances", aggfunc="mean")
juxta_mean.to_csv(out_dir / "misty_juxta_mean.csv")
sns.heatmap(juxta_mean, cmap="Blues", linewidths=0.5)
plt.title(f"MISTy juxta importance, mean of {interactions['sample'].nunique()} samples")
plt.savefig(out_dir / "misty_juxta_mean.png", bbox_inches="tight")
