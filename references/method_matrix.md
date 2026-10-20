# 方法 × 平台适用表

✅ 直接可用　⚠️ 需要调整　❌ 基本不行（原因在最后一列）

| 方法 | 脚本 | 语言 | Visium | HD | Xenium/CosMx | CODEX/IMC | 限制条件 |
|---|---|---|---|---|---|---|---|
| 质控 + 聚类 | preprocess_rna.py / preprocess_protein.py | Py | ✅ | ✅ | ✅ | ✅（蛋白版） | 蛋白不能用 RNA 的标准化 |
| RCTD 去卷积 | export_for_rctd.py → run_rctd.Rmd → import_rctd.py | Py→R→Py | ✅ full | ⚠️ 16µm doublet | ❌ 不需要 | ❌ | 需要单细胞参考 |
| cell2location | run_cell2location.py | Py | ✅ | ⚠️ 太慢 | ❌ | ❌ | 需要参考，最好有 GPU |
| 标签转移 | label_transfer.py | Py | ❌ spot 是混合的 | ✅ 细胞分割后 | ✅ | ⚠️ 需要蛋白与 RNA 对应 | 需要参考 |
| marker 注释 | marker_annotation.py | Py | ⚠️ 只是打分 | ✅ | ✅ | ✅ | marker 要在 panel 里 |
| 空间域 CellCharter | domain_cellcharter.py | Py | ✅ | ✅ | ✅ | ✅ | — |
| 空间域 BANKSY | domain_banksy.Rmd | R | ✅ | ✅ | ✅ | ✅ | 多样本时坐标要错开 |
| niche / RCN / CN | niche_composition.py | Py | ✅ 用细胞比例 | ✅ | ✅ | ✅ | 需要先注释 |
| 邻域富集 / 共现 | nhood_enrichment.py | Py | ⚠️ 只能用主导类型 | ✅ | ✅ | ✅ | — |
| 组间差异邻域 | diff_nhood.py | Py | ⚠️ | ✅ | ✅ | ✅ | 2 组，每组 ≥ 2 个样本（≥ 3 个更可靠） |
| MISTy 共定位 | misty_colocalization.py | Py | ✅ | ✅ | ✅ | ✅（原本为 IMC 开发） | — |
| 空间配受体（LIANA+） | lr_spatial_liana.py | Py | ✅ | ✅ | ⚠️ 先查配受体覆盖 | ❌ | 需要配体和受体都测到 |
| 邻域内配受体（LIANA+ CellPhoneDB） | ligrec_neighborhood.py | Py | ❌ | ✅ | ✅ | ❌ | 需要单细胞级标签；panel 太小会报覆盖率不足 |
| COMMOT 信号流向 | commot_flow.py | Py | ✅ | ⚠️ 内存 | ⚠️ 内存 | ❌ | 单样本运行 |
| 空间 CNV | cnv_infercnvpy.py | Py | ✅ | ⚠️ 用 16µm 或细胞分割 | ❌（5K panel 很勉强） | ❌ | 需要全转录组和基因位置 |
| 空间可变基因 | svg_moran.py | Py | ✅ | ✅ | ✅ | ✅ | — |
| 热点分析 | hotspot_getis.py | Py | ✅ | ✅ | ✅ | ✅ | 任何连续得分都行 |
| NMF 表达程序 | nmf_programs.py | Py | ✅ | ✅（先下采样） | ⚠️ 程序数要少 | ⚠️ | 输入必须非负 |
| 按微环境分亚型 | nmf_subtype.py | Py | ❌ | ✅ | ✅ | ✅ | 需要 nhood_composition |
| 空间趋势 | spatial_tendency_soapy.py | Py | ✅ | ✅ | ✅ | ✅ | 单样本运行 |
| RNA 速度（课程第2课） | 未收录 | Py | ⚠️ 只限 polyA 捕获 | ❌ 探针法 | ❌ | ❌ | 需要剪接 / 未剪接 reads |
| 多组通讯 MultiNicheNet（课程第2课） | 未收录，走 bridge | R | ⚠️ 伪 bulk | ✅ | ⚠️ | ❌ | 多样本多分组 |

总体来看：**基于邻域和空间关系的方法基本通用**（只需要坐标加标签或比例矩阵）；**依赖表达谱的方法受平台限制**（基因数和捕获方式）。
