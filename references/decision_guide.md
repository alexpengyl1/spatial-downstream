# 分析路线推荐：平台 × 有没有单细胞参考

先回答三个问题，就能定下路线：

1. **一个观测点是不是一个细胞？** 是 → `resolution: cell`；不是（spot 里混了多个细胞）→ `resolution: spot`
2. **测了多少基因？** 全转录组 → CNV、通路、NMF 程序都能做；靶向 panel（≤ 5000 个基因）→ 只做基于细胞类型和邻域的分析
3. **有没有配对或同组织的单细胞参考？** 有 → 去卷积或标签转移；没有 → marker 打分 + 人工核对

## 平台归类

| 平台 | 分辨率 | resolution | 基因 | 读入脚本 | 流程 |
|---|---|---|---|---|---|
| 10X Visium（v1/v2、CytAssist） | 55µm spot，1–10 个细胞 | spot | 全转录组 | load_visium.py | P1（有参考）/ P2（无参考） |
| Stereo-seq bin50 / bin100 | 25 / 50µm | spot | 全转录组 | load_stereo.py | P1（有参考）/ P2（无参考） |
| Visium HD 16µm / 8µm bin | 16 / 8µm，1–2 个细胞 | cell | 全转录组 | load_visium_hd.py | P3 / P4 |
| Visium HD 细胞分割、Stereo-seq cellbin、Slide-seq | 近单细胞 | cell | 全转录组（稀疏） | load_visium_hd.py / load_stereo.py | P3 / P4 |
| Xenium、CosMx、MERSCOPE | 单细胞 | cell | 靶向 300–6000 个基因 | load_xenium.py / load_cosmx.py | P5 / P6 |
| CODEX、IMC、mIHC、MIBI | 单细胞 | cell | 蛋白 20–60 个 | load_codex.py | P7 |

## 各场景推荐

### 场景 1（P1）：低精度（Visium）+ 有单细胞参考
| 模块 | 推荐 | 语言 | 备选 |
|---|---|---|---|
| 注释 | RCTD full 模式 | Python → **R** → Python | cell2location（Python，有 GPU 时更好） |
| 空间域 | CellCharter 聚合邻居 | Python | BANKSY（R，domain_banksy.Rmd） |
| niche | 对细胞比例做 KMeans（第5课） | Python | — |
| 共定位 | MISTy（细胞比例） | Python | — |
| 通讯 | LIANA+ bivariate；COMMOT 看信号方向 | Python | SpaGene（R） |
| 个性化 | CNV（infercnvpy）、SVG、热点、NMF 程序、空间趋势 | Python | inferCNV（R，课程第2课） |

### 场景 2（P2）：低精度 + 没有单细胞参考
- 注释：marker 打分（marker_annotation.py），结果是 `obsm['cell_scores']`，**不是比例**，只能当细胞组成的近似。config 里设 `composition_key: cell_scores`。
- 更推荐以"分子 niche"为主：直接用表达聚类（02 的 leiden）和空间域（domain）讲故事，再用 marker 分解释每个域。
- 能做：空间域、MISTy（基于打分）、LIANA+ 配受体、CNV（参考用免疫富集的区域）、SVG、NMF。
- 也可以用无参考去卷积（STdeconvolve，R）得到"主题"比例，再按 marker 命名；本 skill 没收录脚本。

### 场景 3（P3）：Visium HD + 有单细胞参考
- 注释：16µm bin 用 RCTD **doublet** 模式（一个 bin 最多两种细胞），`obs['cell_type']` 取 first_type。
  耗时参考：doublet 模式 1 万个 bin、4 核约 3.5 分钟，几十万个 bin 要按小时算（threads 尽量设大）。run_rctd.Rmd 已为 doublet 模式调低 UMI 阈值。
  8µm 以下的 bin 太稀疏，RCTD 会丢掉很多 bin；最好用 Space Ranger 的细胞分割结果，再走标签转移（label_transfer.py）。
- 之后按单细胞级平台处理：邻域富集、差异邻域、niche、邻域内通讯、空间亚型。
- 数据量大（几十万到上百万个点）：所有脚本都用稀疏邻居图；NMF 先下采样拟合；COMMOT 只取局部区域。

### 场景 4（P4）：Visium HD + 没有参考
- marker_annotation.py 按 Leiden 聚类整簇注释；先看 `markers_leiden.csv`，必要时手动改 `cell_type`。

### 场景 5（P5）：Xenium / CosMx + 有参考
- 注释：标签转移（Harmony + KNN），看 `cell_type_score` 低的细胞是不是集中在某些区域（可能是分割错误或参考里没有的类型）。
- 重点做：RCN（niche_composition）、Delaunay 邻域富集、差异邻域、按微环境分亚型、邻域内配受体。
- 不做：CNV、RNA 速度、COMMOT、通路富集（panel 基因太少）。

### 场景 6（P6）：Xenium / CosMx + 没有参考
- Xenium 自带的 panel 一般有设计好的 marker，marker_annotation.py 按聚类注释即可；其余同场景 5。

### 场景 7（P7）：CODEX / IMC
- 标准化用 arcsinh（preprocess_protein.py），不用 RNA 的流程。
- 注释几乎都是 marker 加人工判断；niche（CN）和邻域分析是核心。
- 不做：去卷积、CNV、配受体通讯。

## 一个完整项目的推荐顺序（所有场景通用）

```
A 读入 → B 质控/聚类 → C 注释 → D 空间域 + niche → E 邻域/共定位（+ 组间差异）→ F 通讯 → G 个性化
01_raw    02_clustered    03_annotated   04_domain → 05_niche   （E/F/G 都读 05_niche.h5ad，结果写成表格和图）
```

每一步只读上一步的 h5ad，所以可以从任意一步重新开始，也可以只跑其中几个模块。
