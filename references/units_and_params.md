# 坐标单位、数据结构约定和关键参数

## 坐标：统一换算成 µm

所有读入脚本都把坐标存成两份：`obsm['spatial']`（平台原始单位，供叠加 HE 图）和 `obsm['spatial_um']`（µm）。
**所有半径、距离参数都按 µm 理解，都用 `spatial_um`。** 课程代码里 `radius=100 / 150 / 500` 没有单位，换平台时容易差出一个数量级。

| 平台 | 原始单位 | 换算 |
|---|---|---|
| Visium | 全分辨率图像像素 | µm/像素 = 55 / `spot_diameter_fullres` |
| Visium HD | 全分辨率图像像素 | `scalefactors_json.json` 里的 `microns_per_pixel` |
| Stereo-seq | DNB 坐标 | × 0.5（DNB 中心距 500nm；按芯片版本核对） |
| Xenium | µm | 不需要换算 |
| CosMx | 像素 | × 0.12028（用 `CenterX/Y_global_px`，不要用 local） |
| CODEX / IMC | 像素 | × 显微镜像素大小（CODEX 常见 0.377；IMC 为 1） |

## 半径怎么定

| 分析 | 低精度（spot） | 单细胞级 |
|---|---|---|
| 直接邻居 | 6 个最近邻（Visium 六边形） | Delaunay，或半径 20–30µm |
| 细胞 niche / RCN / CN | 不需要（spot 组成本身就是邻域） | 半径 30–100µm，或 10–20 个最近邻 |
| 配受体作用范围 | 100–250µm | 30–80µm |
| 共现随距离变化 | 0–1000µm | 0–300µm |

## 统一的数据结构约定

| 位置 | 含义 | 由谁写入 |
|---|---|---|
| `obs['sample']`、`obs['condition']` | 样本、分组 | A_load |
| `layers['counts']` | 原始 counts（R 端和 cell2location 都用它） | A_load |
| `layers['raw']` | 蛋白原始强度 | load_codex.py |
| `X` | log1p 标准化表达（蛋白组为 arcsinh 后的 z-score） | B_preprocess |
| `obsm['X_pca']` / `obsm['X_pca_harmony']` | 降维（多样本用 harmony 版） | B_preprocess |
| `obs['leiden']` | 表达聚类（"分子 niche"） | B_preprocess |
| `obs['cell_type']` | 细胞类型（低精度为主导类型） | C_annotate |
| `obsm['cell_props']` | 细胞类型比例（去卷积） | import_rctd.py / run_cell2location.py |
| `obsm['cell_scores']` | marker 打分 | marker_annotation.py |
| `obs['domain']` | 空间域 | domain_cellcharter.py / domain_banksy.Rmd |
| `obs['niche']`、`obsm['nhood_composition']` | 细胞 niche / RCN / CN 及邻域组成 | niche_composition.py |

文件链：`01_raw → 02_clustered → 03_annotated → 04_domain → 05_niche`，E/F/G 模块都读 `05_niche.h5ad`，只输出表格和图。

## 多样本的两条规则

1. 空间邻居图一律用 `library_key="sample"` 按样本分别构建，不能跨样本连边。
2. 只看坐标的工具（LIANA 的 spatial_neighbors、COMMOT、SOAPy、Getis-Ord）每个样本单独运行。

## 课程代码里的参数，以及本 skill 的改动

| 课程做法 | 问题 | 本 skill |
|---|---|---|
| `sc.pp.regress_out(total_counts, pct_mt)` | 大数据非常慢；和 Harmony 一起用意义不大 | 去掉 |
| 先只保留高变基因，再保存 h5ad | 后面的去卷积、CNV、配受体都缺基因 | 保留全部基因，PCA 时用 `mask_var` |
| Xenium 也选 2000 个高变基因 | 靶向 panel 本来就只有几百个基因 | 基因数 ≤ 1000 时全部使用 |
| `sq.gr.ligrec` 之前先建半径图 | ligrec 不使用空间图；squidpy 1.6.5 上 p 值全为 NaN | 用 LIANA+ bivariate，或先圈出邻域子集再跑 LIANA+ CellPhoneDB |
| `co_occurrence` 的第一个区间直接相减（第16课） | 没有统计检验 | 每个样本做邻域富集，两组之间对 z 分数做 Welch t 检验 |
| `library_key='treatment'`（第23课） | 同组样本之间会连边 | 邻居图一律按样本构建 |
| `scipy.spatial.distance_matrix` 全部两两距离（第14课） | HD 上百万个点会耗尽内存 | 稀疏半径邻居图 × one-hot 矩阵 |
| CosMx 用 `CenterX_local_px` | 不同 FOV 的坐标重叠 | 用 global 坐标 |
