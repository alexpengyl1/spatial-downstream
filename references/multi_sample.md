# 多样本联合，还是只能单样本？

空间分析比单细胞多了一个问题：**哪些步骤可以把所有样本放在一起做，哪些必须一张切片一张切片地做**。判断标准分两层：

> 1. **只用表达、不用坐标的步骤，可以多样本联合**（和单细胞一样，用 Harmony 处理批次）。
> 2. **用到坐标的步骤，"邻居"必须来自同一张切片**；如果统计量还要和"随机"比较（置换检验、Moran's I），最好**逐样本计算，再汇总**。

原因：
- 不同切片的坐标系彼此独立，坐标直接拼在一起，邻居就会跨样本；
- 即使邻居不跨样本，把样本合在一起做置换检验或计算 Moran's I，样本之间细胞组成或表达水平的差异（批次）也会混进结果，而且程序不会报错。

## 四类处理方式

| 类型 | 含义 | 结果怎么跨样本比较 |
|---|---|---|
| **多样本联合** | 所有样本放进一个对象、一个模型；用到坐标时邻居图按样本分别构建 | 所有样本共用同一套聚类、注释、空间域、niche 编号，可以直接比较占比 |
| **逐样本计算，再汇总** | 每个样本单独计算统计量（脚本自动循环），再取平均或统计在几个样本里成立 | 看汇总表；每个样本的结果也单独保存，可按分组比较 |
| **单样本** | 脚本一次只处理一个样本（顶部 `sample = ...`），需要逐个运行 | 自己逐个运行后对比 |
| **需要多样本** | 以样本为重复单位做统计检验 | 本身就是跨样本比较，每组至少 2 个样本 |

## 课程演示代码怎么处理样本，以及本 skill 的结论

下表"演示代码"一列，是从课程原始脚本里实际看到的写法；"结论"一列是本 skill 采用的处理方式。

| 分析（课） | 演示代码怎么处理样本 | 结论 | 原因 |
|---|---|---|---|
| 质控、聚类（1、3、11） | 第1、3课：多个样本 `sc.concat` 合并，`highly_variable_genes(batch_key=...)`，Harmony / BBKNN 去批次，然后在合并数据上聚类；第11课一次读一个 HD 样本 | **多样本联合** | 只用表达；联合聚类后各样本的聚类编号才可比，第3课课件也强调"批次校正，确保不同样本的细胞类型可比性" |
| 去卷积 RCTD（4、5、13） | 每次运行读一个空间样本（`--name`）；课件写"RCTD 支持单样本分析" | **逐样本计算** | RCTD 每次拟合一张切片的平台效应；各样本共用同一份单细胞参考，所以细胞类型一致，结果可以直接合并 |
| 去卷积 cell2location（4） | 参考数据设 `batch_key='Sample'`；课件写"支持多样本整合" | **多样本联合** | 模型里把样本作为批次一起训练 |
| 标签转移（20） | 单细胞和 Xenium 合并后，用 Harmony 按 `dataset` 去批次，再做 KNN | **多样本联合** | 只用表达 |
| 细胞 niche（5） | 所有 spot 的 RCTD 比例一起做 SNN 聚类；再按病人（`orig.ident`）取各 niche 细胞比例的中位数，用 Wilcoxon 检验每个 niche 的特征细胞 | **多样本联合**，并输出按样本的组成表 | niche 要跨样本可比，所以联合聚类；但 niche 的特征细胞要以病人 / 样本为单位检验，所以本 skill 另外输出 `niche_composition_by_sample.csv` |
| 共定位 MISTy（6） | 演示只读一个 `Visium.h5ad`；课件标题"单样本多样本共定位"，并以"共定位分析（多样本比较）"引用心肌梗死空间多组学文章 | **逐样本计算，再汇总** | MISTy 只看坐标；那篇文章就是每张切片单独建模，再比较各组的重要性 |
| 通讯 COMMOT、SpaGene（7） | 一次读一个样本（`read_visium(library_id=sample)`、一个 rds） | **单样本** | 最优传输 / 配受体共表达都在一张切片的距离上计算 |
| 空间 CNV（8） | 读一个 h5ad，用 `regions` 列指定正常参考 | **多样本联合** | 只用表达和基因位置；所有样本共用一组正常参考，CNV 聚类才可比 |
| 热点 spottedpy（9） | `obs['batch'] = obs['sample']`，`create_hotspots(..., relative_to_batch=True)` | **逐样本计算** | 演示代码明确按样本（batch）分别计算热点 |
| 热点 Hotspot 模块分（9） | 读合并好的 `combined_harmony.h5ad`，在全部 spot 的坐标上建 KNN | 演示代码有问题 → **逐样本** | 多样本合并后直接按坐标找邻居会跨样本；本 skill 的 `hotspot_getis.py` 逐样本计算 |
| 空间趋势 SOAPy（10） | 读一个 `test.h5ad`，按 niche 生成掩膜 | **单样本** | 掩膜画在一张切片上 |
| 空间域 BANKSY（12） | 课题目为"依据邻域识别空间域与多样本整合"，演示代码读一个样本 | **多样本联合** | 本 skill 按 BANKSY 官方多样本做法错开坐标后联合计算；CellCharter 按样本建邻居图后联合聚类 |
| niche（14，scNiche） | 在全部细胞上算两两距离；富集按 `library_key='Group'` | **多样本联合**（邻居不跨样本） | 演示代码的两两距离在多样本时会跨样本，本 skill 改为按样本建邻居图 |
| 邻域富集（15、19、22） | 第15、22课：在一个对象上建邻居图、做富集，没有样本参数（按单样本写）；第19课：`for sample in adata.obs.fov.unique()` 逐个 FOV 做 `nhood_enrichment` | **逐样本计算，再汇总** | 置换检验打乱的是细胞类型标签，多个样本合在一起打乱会把样本间组成差异混进零分布；第21课课件也说"一般邻域分析侧重单样本内的空间结构" |
| 邻域差异（16、23） | 第16课：两个样本各自算 `co_occurrence` 再相减；第23课：按 `treatment` 分组，`library_key='treatment'` | **需要多样本** | 以样本为重复单位；第23课把整组当一个 library 会让同组样本之间连边，本 skill 改为按样本 |
| 配受体（17、25） | 一个对象上做 `ligrec`；第25课在一个对象里按半径找邻居 | 空间配受体 **逐样本**；邻域内配受体 **联合** | LIANA+ 的局部得分依赖切片内的空间核；`ligrec` 的置换只打乱细胞类型、不看坐标，邻域圈好后可以合并 |
| NMF（18、24） | 第18课：一个 h5ad 上做 NMF；第24课：在 scimap 邻域计数上做 NMF，再 `pie(group_by='treatment')` 按分组比较 | **多样本联合** | 只用表达或已经按样本算好的邻域组成；联合后程序 / 亚型编号可比 |
| RCN（21） | `spatial_count(..., imageid='fov')` 按 FOV 统计邻域，再一起聚类；`pie(group_by='fov')`；课件："RCN 关注跨样本 / 区域的复现性模式，需要多个样本或不同组织区域" | **多样本联合**（最好有多个样本） | 邻域按样本统计，聚类联合做，才能找出在多个样本里反复出现的 RCN |
| CN（27，SOAPy） | `make_network(sample_key='fov')`、`get_c_niche(sample_key='fov')` | **多样本联合**（邻居不跨样本） | 网络按 FOV 分别建，CN 联合聚类；同课 `pro.RCN.py` 在全部细胞上按半径找邻居，多 FOV 时会跨 FOV |
| 空间可变基因（本 skill 新增，替代第7课 SpaGene） | — | **逐样本计算，再汇总** | 实测：一个没有空间模式、只在 Treated 组整体升高的基因，合并计算 Moran's I = 0.19（会被当成空间可变基因），逐样本计算都约为 0 |

## 本 skill 每个脚本的处理方式

| 脚本 | 类型 | 输出 |
|---|---|---|
| `A_load/*` | 逐样本读入 → 合并 | `01_raw.h5ad`（条形码加样本后缀） |
| `preprocess_rna.py` / `preprocess_protein.py` | 多样本联合（Harmony） | `02_clustered.h5ad` |
| `run_rctd.Rmd` | 逐样本 | `rctd/<样本>/weights.csv` → `03_annotated.h5ad` |
| `run_cell2location.py`、`label_transfer.py`、`marker_annotation.py` | 多样本联合 | `03_annotated.h5ad` |
| `domain_cellcharter.py` / `domain_banksy.Rmd` | 多样本联合（邻居不跨样本） | `04_domain.h5ad` |
| `niche_composition.py` | 多样本联合（邻居不跨样本） | `05_niche.h5ad`、`niche_by_sample.csv`、`niche_composition_by_sample.csv` |
| `nhood_enrichment.py` | 逐样本 → 平均 | `nhood_enrichment_<样本>.csv`、`nhood_enrichment_mean.csv/png`（格子里标 z > 0 的样本数） |
| `diff_nhood.py` | 需要多样本（≥ 2 组 × ≥ 2 样本） | `diff_nhood_zscore.csv`、`diff_nhood_pvalue.csv` |
| `misty_colocalization.py` | 逐样本 → 平均 | `misty_interactions_by_sample.csv`、`misty_juxta_mean.csv/png` |
| `lr_spatial_liana.py` | 逐样本 | `lr_global_<样本>.csv` |
| `ligrec_neighborhood.py` | 多样本联合（邻居不跨样本） | `ligrec_<细胞类型>.csv` |
| `commot_flow.py` | 单样本 | `commot_<通路>.png` |
| `cnv_infercnvpy.py` | 多样本联合 | `cnv_per_spot.csv` |
| `svg_moran.py` | 逐样本 → 汇总 | `svg_moran_<样本>.csv`、`svg_moran.csv`（平均 I、在几个样本里显著） |
| `nmf_programs.py`、`nmf_subtype.py` | 多样本联合 | `nmf_top_genes.csv`、`subtype_markers.csv`、`subtype_by_condition.csv` |
| `hotspot_getis.py` | 逐样本 | `hotspot_distance.csv` |
| `spatial_tendency_soapy.py` | 单样本 | `tendency_*.csv` |

## 实际怎么用

1. **把所有样本写进同一个 `samples.tsv`，一次跑完整套流程。** 联合步骤会自动跨样本，逐样本步骤会自动循环。不要每个样本单独建一个项目，否则聚类、niche、空间域的编号在样本之间对不上。
2. **看逐样本步骤的结果时，先看汇总，再看一致性。** 比如邻域富集热图里每个格子标了"z > 0 的样本数 / 样本总数"，只在一个样本里出现的关系要谨慎解释。
3. **比较分组（Control vs Treated）**：
   - 组成差异：`niche_by_sample.csv`、`niche_composition_by_sample.csv`、`subtype_by_condition.csv`，以样本为单位做检验；
   - 邻域关系差异：`diff_nhood.py`，每组至少 2 个样本；
   - 共定位差异：`misty_interactions_by_sample.csv`，按分组比较各样本的重要性。
4. **只有 1 个样本**：所有"联合"和"逐样本"步骤照常运行（汇总就等于这个样本本身）；`diff_nhood.py` 跳过（流程里设 `HAS_GROUPS=no`）。
5. **单样本步骤（COMMOT、空间趋势）**：先对代表性样本运行，需要时再改 `sample` 逐个运行。
6. **自己新加步骤时**，按顺序问自己：
   - 这一步用坐标吗？不用 → 可以联合；
   - 用的话，邻居会不会跨样本？会 → 按样本建邻居图，或在样本上循环；
   - 结果要和随机比较（置换检验、自相关）吗？要 → 逐样本计算再汇总。
