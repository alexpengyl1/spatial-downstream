---
name: spatial-downstream
description: 空间转录组 / 空间蛋白组下游分析流程（Visium、Visium HD、Stereo-seq、Slide-seq、Xenium、CosMx、MERSCOPE、CODEX、IMC）。按平台精度（低精度 spot / 高精度单细胞级 / 蛋白组）和有无单细胞参考推荐路线，并给出每一步用 Python 还是 R、R 与 Python 之间怎么交换数据。覆盖读入、质控聚类、去卷积或注释（RCTD、cell2location 等）、空间域（CellCharter、BANKSY）、niche / RCN / CN、邻域富集与组间差异、MISTy 共定位、空间配受体通讯（LIANA+、squidpy、COMMOT）、空间 CNV、SVG、热点、NMF、空间亚型、空间趋势。Use for spatial transcriptomics downstream analysis, deconvolution, niches, neighbourhoods, spatial cell-cell communication, or R + Python spatial pipelines.
---

# 空间转录组下游分析

原始示例脚本与本 skill 脚本的对应关系，以及对原始代码做的修正，见 `references/script_map.md`。
所有脚本按 readable-coding 规范写：从上到下运行，用 `# %%` 分块（可在 Jupyter / VS Code 里逐块运行），每块开头有"步骤 / 输入 / 输出"三行注释，参数集中在顶部或 `config.yaml`。

## 第一步：确定路线

问清楚（或从数据判断）三件事，再查 `references/decision_guide.md`：

1. **平台和分辨率**：spot 级（Visium、Stereo-seq bin≥50）→ `resolution: spot`；单细胞级（HD、cellbin、Xenium、CosMx、CODEX）→ `resolution: cell`
2. **基因数**：全转录组，还是靶向 panel / 蛋白 panel（决定 CNV、通路、配受体能不能做，见 `references/method_matrix.md`）
3. **有没有单细胞参考**：决定注释方法

| 编号 | 场景 | 平台 | 注释方法 | R 步骤 |
|---|---|---|---|---|
| P1 `pipelines/P1_lowres_scrna.sh` | 低精度 + 单细胞参考 | Visium、Stereo-seq bin≥50 | RCTD full（或 cell2location） | run_rctd.Rmd |
| P2 `pipelines/P2_lowres_noref.sh` | 低精度，无参考 | 同上 | marker 打分（cell_scores） | 无 |
| P3 `pipelines/P3_hd_scrna.sh` | 高精度测序型 + 参考 | Visium HD、Stereo-seq bin≤20 / cellbin、Slide-seq | RCTD doublet | run_rctd.Rmd |
| P4 `pipelines/P4_hd_noref.sh` | 高精度测序型，无参考 | 同上 | 按聚类 marker 注释 | 无 |
| P5 `pipelines/P5_imaging_scrna.sh` | 成像型单细胞 + 参考 | Xenium、CosMx、MERSCOPE | Harmony + KNN 标签转移 | 无 |
| P6 `pipelines/P6_imaging_noref.sh` | 成像型单细胞，无参考 | 同上 | 按聚类 marker 注释 | 无 |
| P7 `pipelines/P7_protein.sh` | 空间蛋白组 | CODEX、IMC、mIHC、MIBI | marker + 人工核对 | 无 |

每套流程跑哪些步骤，见 `README.zh-CN.md` 第 3 节的对照表；流程图在 `docs/`，分三层：图 1 `decision_zh.svg`（根据数据选出 P1–P7 中的一套）→ 图 2 `steps_zh.svg`（七套流程每一步用什么、哪些不适用）→ 图 3 `P1…P7_zh.svg`（单套流程的详细步骤，每一步标出读什么、写什么）；图 4 `dataflow_zh.svg` 说明 01–05 各阶段文件的字段和下游脚本读哪一部分；英文版为 `_en`；向用户介绍流程时可以直接给图。改了流程后用 `python docs/make_figures.py` 重新生成。

## 第二步：建项目

```bash
SKILL_DIR=~/.claude/skills/spatial-downstream      # 在 Codex 里是 ~/.codex/skills/spatial-downstream
cp "$SKILL_DIR"/assets/configs/P1_lowres_scrna.yaml config.yaml    # 换成对应流程的预设
cp "$SKILL_DIR"/assets/samples_template.tsv samples.tsv
cp "$SKILL_DIR"/assets/markers_template.yaml markers.yaml
bash "$SKILL_DIR"/pipelines/P1_lowres_scrna.sh
```

- 改 `config.yaml`：species、scrna_reference、radius_um、n_niches 等；字段定义见 `schemas/config.schema.json`，每套流程第一步会用它检查配置
- 改 `samples.tsv`：每行一个样本（sample、condition、path）
- 流程脚本顶部的变量（`PY`、`R_CMD`、`LOADER`、`HAS_GROUPS`）可以改，也可以用环境变量覆盖
- 环境：`envs/spatial-py.yml`（固定版本，已测试）、`envs/install_r_packages.R`（RCTD / BANKSY，装到系统 R；bioconda 上两者的 R 版本对不上）；cell2location、COMMOT、stereopy、SOAPy 各用单独环境（见 yml 末尾）
- 每次运行最后会写 `results/provenance/`（配置副本、软件版本、输入输出校验值），交付结果时一起给用户

## 第三步：按模块运行

所有脚本都在**项目目录**下运行（读 `./config.yaml`）。每个模块读上一个模块的 h5ad，所以可以从任意一步重跑。

| 模块 | 脚本（`scripts/`） | 语言 | 输入 → 输出 |
|---|---|---|---|
| A 读入 | `A_load/load_{visium,visium_hd,stereo,xenium,cosmx,codex}.py` | Py | 原始输出 → `01_raw.h5ad` |
| B 质控聚类 | `B_preprocess/preprocess_rna.py`；蛋白用 `preprocess_protein.py` | Py | → `02_clustered.h5ad` |
| C 注释 | `C_annotate/prepare_reference.py`（有参考时先运行）<br>`export_for_rctd.py` → `run_rctd.Rmd` → `import_rctd.py`<br>或 `run_cell2location.py` / `label_transfer.py` / `marker_annotation.py` | Py / **R** | → `03_annotated.h5ad` |
| D 空间域 + niche | `D_domain_niche/domain_cellcharter.py`（或 `bridge/h5ad_to_mtx.py` → `domain_banksy.Rmd` → `bridge/import_obs_csv.py`）<br>`niche_composition.py` | Py / **R** | → `04_domain.h5ad` → `05_niche.h5ad` |
| E 邻域 / 共定位 | `E_neighborhood/nhood_enrichment.py`（单细胞级）<br>`diff_nhood.py`（组间差异）<br>`misty_colocalization.py`（低精度首选） | Py | `05_niche.h5ad` → 图和表 |
| F 通讯 | `F_communication/lr_spatial_liana.py`（所有 RNA 平台）<br>`ligrec_neighborhood.py`（单细胞级）<br>`commot_flow.py`（低精度、看方向） | Py | `05_niche.h5ad` → 图和表 |
| G 个性化 | `G_special/cnv_infercnvpy.py`、`svg_moran.py`、`hotspot_getis.py`、`nmf_programs.py`、`nmf_subtype.py`、`spatial_tendency_soapy.py` | Py | `05_niche.h5ad` → 图和表 |

## 使用示例：4 个 Visium 样本 + 单细胞参考

用户的请求：

> 我有 4 个人肿瘤的 Visium 样本，2 个对照、2 个处理，还有一份同组织的单细胞数据 `ref/scrna.h5ad`（原始 counts，细胞类型在 `cell_type` 列）。想看处理之后肿瘤周围的免疫细胞组成有没有变化。

**1. 选路线。** Visium 是 spot 级数据，又有单细胞参考，所以用 P1：RCTD 去卷积得到每个 spot 的细胞比例，再按比例找 niche。回复时用一句话说明为什么不是 P2（有参考就不用 marker 打分）。

**2. 建项目。**

```bash
SKILL_DIR=~/.claude/skills/spatial-downstream      # 在 Codex 里是 ~/.codex/skills/spatial-downstream
mkdir tumor_visium && cd tumor_visium
cp "$SKILL_DIR"/assets/configs/P1_lowres_scrna.yaml config.yaml
cp "$SKILL_DIR"/assets/samples_template.tsv samples.tsv
cp "$SKILL_DIR"/assets/markers_template.yaml markers.yaml
```

`samples.tsv` 每行一个样本，`path` 指向 Space Ranger 的 `outs/` 目录：

```
sample	condition	path
C1	Control	/data/visium/C1/outs
C2	Control	/data/visium/C2/outs
T1	Treated	/data/visium/T1/outs
T2	Treated	/data/visium/T2/outs
```

`config.yaml` 只需核对 `species: human`、`scrna_reference: ref/scrna.h5ad`、`scrna_celltype_key: cell_type`；`n_niches` 先用默认的 8，看过 niche 图再调。

**3. 运行。** `bash "$SKILL_DIR"/pipelines/P1_lowres_scrna.sh`。跑完先看 `results/qc_violin.png`，过滤阈值不合适就改 `preprocess_rna.py` 顶部的参数，从这一步重跑。

**4. 用哪些结果回答问题。**

| 问题 | 看哪里（都在 `results/`） |
|---|---|
| 每个 spot 里有哪些细胞 | `03_annotated.h5ad` 的 `obsm['cell_props']`；RCTD 报告 `reports/run_rctd.html` |
| 有哪几类微环境，各由什么细胞组成 | `niche_composition_heatmap.png`、`spatial_niche.png` |
| 各样本里每个 niche 占多少，对照和处理比 | `niche_by_sample.csv`、`niche_by_sample.png` |
| 同一个 niche 在不同样本里的细胞组成 | `niche_composition_by_sample.csv` |
| 哪些细胞类型挨在一起 | `misty_juxta_mean.png`（各样本的平均） |
| 肿瘤和免疫细胞之间的配受体 | `lr_local_top.png`、`lr_global_<样本>.csv` |

**5. 向用户说明的几点。**
- niche 是 4 个样本一起聚类的，niche 编号在样本之间可以直接比较；RCTD、MISTy 和 Moran's I 是逐样本算的，表里给的是平均值。
- 每组只有 2 个样本，`niche_by_sample.csv` 只能看趋势。要做统计检验，每组至少 3 个样本更可靠。
- spot 不是单个细胞，所以 P1 不跑邻域富集和组间邻域差异（`nhood_enrichment.py`、`diff_nhood.py`），共定位看 MISTy。
- 交付时把 `results/provenance/` 一起给用户。

## 笔记本和报告

- 每个 `.py` 在 `notebooks/` 下有同名 JupyterLab 笔记本（代码相同，加了小节标题和结果说明）；改了 `.py` 后运行 `python tools/make_notebooks.py` 同步。用户想逐格运行时给笔记本，流程化运行用 `.py`。
- R 步骤是 R Markdown 报告（`.Rmd`）：在 RStudio 里 Knit，或由流程用 `rmarkdown::render()` 渲染到 `results/reports/`。注意 `render()` 的参数是惰性求值的，先 `project <- getwd()` 再传 `params = list(project_dir = project)`。

## R 与 Python 联合（详见 `references/r_python_bridge.md`）

- 主线对象只有 h5ad；R 端只读交换文件、只返回结果表（CSV，第一列为条形码）
- 交换格式：`counts.mtx`（基因 x 细胞，原始整数 counts）+ `genes.txt` + `barcodes.txt` + `obs.csv` + `coords.csv`（µm）
- 读回时一律按条形码对齐（`.reindex(adata.obs_names)`），R 端丢掉的细胞记为 NaN
- 现成的桥接脚本：`bridge/h5ad_to_mtx.py`、`bridge/import_obs_csv.py`、`bridge/seurat_to_mtx.Rmd`、`bridge/mtx_to_h5ad.py`
- 交互式分析时可以用 `zellkonverter::readH5AD()`；流程化运行时仍用文件交换

## 必须遵守的约定（详见 `references/units_and_params.md` 和 `schemas/anndata_contract.md`）

1. **坐标单位**：所有半径、距离都用 `obsm['spatial_um']`（µm）。各平台的换算在读入脚本里完成。
2. **多样本**：只用表达的步骤多样本联合；用坐标的步骤邻居一律不跨样本；要和随机比较的统计量（邻域富集、MISTy、Moran's I、LIANA 空间核、Getis-Ord、共现、RCTD）逐样本计算再汇总；COMMOT 和 SOAPy 只做单样本；组间差异需要每组 ≥ 2 个样本。每一步属于哪一类、课程演示代码怎么写、为什么，见 `references/multi_sample.md`，向用户说明流程时要讲清楚。
3. **counts**：`layers['counts']` 始终保留原始 counts；`X` 是 log 标准化值；保存的 h5ad 保留全部基因（PCA 只用高变基因）。
4. **低精度平台不要用 nhood_enrichment 当主结果**：spot 不是细胞，用细胞比例做 niche 和 MISTy。
5. **`sq.gr.ligrec` 不使用空间坐标**：要空间约束，就用 `lr_spatial_liana.py`，或先用 `ligrec_neighborhood.py` 圈出邻域。
6. **靶向 panel 不做** CNV、通路富集、RNA 速度；配受体分析之前先查 panel 覆盖了多少配受体对。
7. **多进程**：在 macOS / Windows 上用 `python 脚本.py` 运行时，Python 会用 spawn 方式开子进程，子进程会把整个脚本重新执行一遍，导致报错（EOFError / "bootstrapping phase"）或卡住。
   所以脚本里 squidpy 的置换检验都设了 `show_progress_bar=False`；infercnvpy 一定会开进程池，所以 cnv_infercnvpy.py 开头用 `multiprocessing.set_start_method("fork", force=True)`。
   新加的步骤如果用到进程池（ProcessPoolExecutor、tqdm 的 process_map 等），也按这两种方式处理；在 Jupyter 里运行不受影响。

## 测试

改了脚本之后运行 `bash tests/run_tests.sh <工作目录>`：用各平台真实格式的模拟数据把 P1–P7 全部跑一遍，打印 PASS / FAIL。

## 给用户交付代码时

- 按 readable-coding 规范：注释用用户的语言；只给用户要的步骤；为数据另外选的参数（阈值、半径、分辨率、niche 数）在回复里说明原因和修改位置。
- 先看数据再定阈值（`qc_violin.png`），不要直接套默认值。
- 低质量区域、分割错误、参考里没有的细胞类型，会让注释置信度（`cell_type_score`）偏低，要提醒用户检查。
