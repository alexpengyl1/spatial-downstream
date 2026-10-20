# spatial-downstream

空间转录组和空间蛋白组的下游分析流程。每一步都写明用什么方法、用 Python 还是 R，以及能不能把多个样本放在一起分析。

[English README](README.md) · 版本 1.0.0 · 许可证：MIT · 语言：Python + R

关键词：空间转录组、空间蛋白组、Visium、Visium HD、Stereo-seq、Slide-seq、Xenium、CosMx、MERSCOPE、CODEX、IMC、去卷积、RCTD、cell2location、标签转移、空间域、BANKSY、CellCharter、细胞 niche、细胞邻域、邻域富集、MISTy、配体-受体、LIANA+、COMMOT、空间 CNV、scanpy、squidpy、FAIR。

## 1. 为什么做这个 skill

做空间组学下游分析，最常卡住的是两件事。

一是工具太多，还分散在 R 和 Python 两边。同一个任务往往有好几个候选：去卷积用 RCTD 还是 cell2location？空间域用 BANKSY 还是 CellCharter？通讯用 squidpy、LIANA+、COMMOT 还是 CellChat？有的只有 R 版，有的只有 Python 版。哪个适合自己的平台、两种语言之间怎么传数据，往往说不清楚。

二是不知道哪些分析能多样本一起做。有的步骤可以把所有样本放进一个模型，有的只在一张切片内部才有意义。搞混了程序也不报错，只是"邻居"会跨到别的切片上，结果就错了。

这个 skill 按数据类型和有没有单细胞参考，分成 P1 到 P7 七套流程，用图 1 回答三个问题就能选出一套。每一步只给一个推荐方法，并写明用 Python 还是 R；需要 R 的步骤配有测试过的文件交换方式，数据传过去再传回来。图 2 到图 4 标出每一步读上一步的哪部分数据、写出什么。每一步还标明属于"多样本联合 / 逐样本计算再汇总 / 只能单样本 / 需要多个样本"中的哪一类，并说明原因（见下文"多样本还是单样本？"）。脚本都很短，环境版本固定，每次运行都会留下溯源记录；七套流程都用各平台真实文件格式的模拟数据从头到尾跑通过。

每一步都是从上到下运行的短脚本，用 `# %%` 分块，可以在 Jupyter / VS Code 里逐块运行，也可以 `python 脚本.py` 直接运行。步骤之间按统一的 AnnData 数据约定传递数据；`SKILL.md` 告诉 Claude Code 或 Codex 助手怎么选流程、怎么按数据调整。

## 2. 选流程

文档里有四张图，从粗到细：

| 图 | 文件 | 回答什么问题 |
|---|---|---|
| 图 1 · 选择流程 | `docs/decision_zh.svg` | 我的数据用 P1–P7 中的哪一套？ |
| 图 2 · 分步总览 | `docs/steps_zh.svg` | 每种数据每一步用什么？（七套流程放在一张图里） |
| 图 3 · 单套流程详图 | `docs/P1_zh.svg` … `docs/P7_zh.svg` | 选定的这一套具体怎么跑、按什么顺序运行哪些脚本？每一步读什么、写什么？ |
| 图 4 · 数据流 | `docs/dataflow_zh.svg` | 每一步用上一步文件里的哪一部分、新增哪些字段？（所有流程通用） |

先用图 1 选出流程，在图 2 里找到它那一列，再打开它的图 3 照着运行。图 4 说明串起各步骤的文件里有什么。

**图 1 · 选择流程**（从上往下依次回答三个问题）：

![图 1：选择流程](docs/decision_zh.svg)

| 编号 | 数据类型 | 平台 | 单细胞参考 | 注释方法 | 语言 |
|---|---|---|---|---|---|
| **P1** `P1_lowres_scrna` | spot 级（一个点含多个细胞） | Visium、Stereo-seq bin ≥ 50 | 有 | RCTD full（或 cell2location） | Python + R（RCTD 这一步） |
| **P2** `P2_lowres_noref` | spot 级 | Visium、Stereo-seq bin ≥ 50 | 无 | marker 打分 | Python |
| **P3** `P3_hd_scrna` | 近单细胞，测序型 | Visium HD 8/16 µm、Stereo-seq bin ≤ 20 / cellbin、Slide-seq | 有 | RCTD doublet | Python + R（RCTD 这一步） |
| **P4** `P4_hd_noref` | 近单细胞，测序型 | 同 P3 | 无 | 按聚类 marker 打分 | Python |
| **P5** `P5_imaging_scrna` | 单细胞，靶向基因 panel | Xenium、CosMx、MERSCOPE | 有 | Harmony + KNN 标签转移 | Python |
| **P6** `P6_imaging_noref` | 单细胞，靶向基因 panel | 同 P5 | 无 | 按聚类 marker 打分 | Python |
| **P7** `P7_protein` | 单细胞，蛋白 panel | CODEX/PhenoCycler、IMC、mIHC、MIBI | 不适用 | marker 打分 + 人工核对 | Python |

三个问题依次是：
1. 一个观测点是不是一个细胞？不是 → P1/P2；是 → P3–P7。
2. 是测序型全转录组（P3/P4）、靶向 RNA panel（P5/P6），还是蛋白（P7）？
3. 有没有配对或同组织、带原始 counts 的单细胞参考？有 → 奇数编号；没有 → 偶数编号。

各种情况的细节和注意事项见 [`references/decision_guide.md`](references/decision_guide.md)。

## 3. 每套流程跑哪些步骤

**图 2 · 七套流程的分步总览。** 每一列是一种数据，每一行是一个步骤，按运行顺序从上往下排。实线框默认运行，虚线框是可选步骤，灰色框表示不适用（框里写了原因），红色"R"标记表示 R 步骤。

![图 2：分步总览](docs/steps_zh.svg)

同样的内容用表格列出：

| 模块 | 步骤（脚本） | P1 | P2 | P3 | P4 | P5 | P6 | P7 | 语言 |
|---|---|---|---|---|---|---|---|---|---|
| Z | `validate_config.py` 检查配置 | ● | ● | ● | ● | ● | ● | ● | Py |
| A | `load_visium` / `load_visium_hd` / `load_stereo` / `load_xenium` / `load_cosmx` / `load_codex` | ● | ● | ● | ● | ● | ● | ● | Py |
| B | `preprocess_rna.py`（P7 用 `preprocess_protein.py`） | ● | ● | ● | ● | ● | ● | ● | Py |
| C | `prepare_reference` → `export_for_rctd` → **`run_rctd.Rmd`** → `import_rctd` | ● | | ● | | | | | Py → **R** → Py |
| C | `label_transfer.py` | | | | | ● | | | Py |
| C | `marker_annotation.py` | | ● | | ● | | ● | ● | Py |
| D | `domain_cellcharter.py`（或 **`domain_banksy.Rmd`**） | ● | ● | ● | ● | ● | ● | ● | Py（或 R） |
| D | `niche_composition.py`（niche / RCN / CN） | ● | ● | ● | ● | ● | ● | ● | Py |
| E | `misty_colocalization.py` | ● | ● | | | | | ● | Py |
| E | `nhood_enrichment.py` | | | ● | ● | ● | ● | ● | Py |
| E | `diff_nhood.py`（≥ 2 组、每组 ≥ 2 个样本） | | | ● | ● | ● | ● | ● | Py |
| F | `lr_spatial_liana.py` | ● | ● | ● | ● | ● | ● | | Py |
| F | `ligrec_neighborhood.py`、`commot_flow.py` | 可选 | 可选 | 可选 | 可选 | 可选 | 可选 | | Py |
| G | `svg_moran.py` | ● | ● | ● | ● | ● | ● | ● | Py |
| G | `nmf_programs.py` | ● | ● | ● | ● | | | | Py |
| G | `cnv_infercnvpy.py`、`hotspot_getis.py`、`nmf_subtype.py`、`spatial_tendency_soapy.py` | 可选 | 可选 | 可选 | 可选 | 可选 | 可选 | 可选 | Py |
| Z | `record_provenance.py` 记录溯源 | ● | ● | ● | ● | ● | ● | ● | Py |

● 默认运行；可选 = 在流程脚本里被注释掉，设好参数后取消注释即可。

### 图 3 · 单套流程详图

下面每张图是图 2 中的一列画成的详图。每个步骤下方，"读"是它用到的上一步文件里的哪一部分，"写"是它写出的文件。实心编号的步骤默认运行，虚线编号的是可选步骤。蓝色是 Python 步骤，红色是 R 步骤，黄色是关键输出和单细胞参考，灰色用于配置、溯源和不适用的步骤。重新生成这些图：`python docs/make_figures.py`。

<details><summary><b>P1</b> · 低精度 + 单细胞参考</summary>

![P1](docs/P1_zh.svg)

</details>

<details><summary><b>P2</b> · 低精度，无参考</summary>

![P2](docs/P2_zh.svg)

</details>

<details><summary><b>P3</b> · 高精度测序型 + 参考</summary>

![P3](docs/P3_zh.svg)

</details>

<details><summary><b>P4</b> · 高精度测序型，无参考</summary>

![P4](docs/P4_zh.svg)

</details>

<details><summary><b>P5</b> · 成像型单细胞 + 参考</summary>

![P5](docs/P5_zh.svg)

</details>

<details><summary><b>P6</b> · 成像型单细胞，无参考</summary>

![P6](docs/P6_zh.svg)

</details>

<details><summary><b>P7</b> · 空间蛋白组</summary>

![P7](docs/P7_zh.svg)

</details>

### 图 4 · 数据流

每一步读上一步 h5ad 文件里的某一部分，写出下一个文件，下游脚本都读 `05_niche.h5ad`。字段定义见 [`schemas/anndata_contract.md`](schemas/anndata_contract.md)。

![图 4：数据流](docs/dataflow_zh.svg)

### 多样本还是单样本？

下面的分类对照原始示例代码处理样本的方式核对过。只用表达量的步骤可以多样本联合，批次由 Harmony 校正。用到坐标的步骤，"邻居"必须来自同一张切片。要和随机情况比较的统计量（置换检验、Moran's I）逐样本计算再汇总，因为合在一起算会把样本间的差异混进去。实际使用时，把所有样本写进同一个 `samples.tsv`，整套流程跑一次即可。每个原始示例脚本怎么处理样本、为什么这样分类，见 [`references/multi_sample.md`](references/multi_sample.md)。

| 处理方式 | 步骤 | 原因 |
|---|---|---|
| 多样本联合（一个模型，聚类、niche、空间域编号跨样本可比） | 质控与聚类（Harmony）、cell2location、标签转移、marker 注释、CNV、空间域（CellCharter / BANKSY）、niche / RCN / CN、邻域内配受体、NMF 程序、微环境亚型 | 只用表达量的步骤直接合并。用坐标的步骤先按样本建邻居图，再一起聚类，这样 niche 和空间域的编号在各样本里含义相同。原始的 RCN 示例也是这样做的：按 FOV 统计邻域，再一起聚类。 |
| 逐样本计算，再汇总（脚本自动循环，各样本结果也单独保存） | RCTD、邻域富集（平均 z 值，以及几个样本方向一致）、MISTy（平均重要性）、空间可变基因（平均 Moran's I，以及在几个样本里显著）、共现、LIANA+ 空间配受体、热点与距离 | 统计量只在一张切片内有意义，或者合在一起算会把样本间差异混进去。测试中，一个只有批次差异、没有空间模式的基因，合并计算 Moran's I = 0.19，逐样本计算约为 0。 |
| 单样本（改脚本顶部的 `sample`，逐个运行） | COMMOT 信号流向、空间趋势（SOAPy） | 最优传输和区域掩膜只在一张切片内有意义。 |
| 需要多个样本 | 组间邻域差异（`diff_nhood.py`）；用 `niche_composition_by_sample.csv` 做样本水平的检验 | 样本是重复单位，至少要 2 组、每组 2 个样本；每组 3 个以上更可靠。 |

## 4. 快速开始

### 安装 skill（Claude Code 和 Codex 都可以）

这个 skill 是一个顶层带 `SKILL.md` 的文件夹，Claude Code 和 OpenAI Codex 都读这种格式。从 GitHub 克隆到所用工具的 skills 目录（需要 `git`）：

```bash
# Claude Code（所有项目可用）
git clone https://github.com/alexpengyl1/spatial-downstream.git ~/.claude/skills/spatial-downstream

# Codex
git clone https://github.com/alexpengyl1/spatial-downstream.git ~/.codex/skills/spatial-downstream
```

只想在某个项目里用 Claude Code，就克隆到 `<项目>/.claude/skills/spatial-downstream`。没有 `git` 的话，在 GitHub 页面点 Code → Download ZIP，解压后把文件夹改名为 `spatial-downstream`，放进上面的某个 skills 目录。以后更新运行 `git -C ~/.claude/skills/spatial-downstream pull`。

两个工具都用的话，装一份，再链接到另一个，例如 `ln -s ~/.claude/skills/spatial-downstream ~/.codex/skills/spatial-downstream`。装好后重启工具。它通过 `SKILL.md` 里的 `name` 和 `description` 找到这个 skill，之后直接说"用单细胞参考分析我的 Visium 数据"之类就行。

流程脚本根据自身位置找到其他脚本，装在哪个目录都能运行。只有在流程之外单独运行 `validate_config` 或 `record_provenance`（脚本或笔记本），并且装的位置不是 `~/.claude/skills/spatial-downstream` 时，才需要把 `SPATIAL_SKILL_DIR` 设成安装目录。

### 运行流程

```bash
SKILL_DIR=~/.claude/skills/spatial-downstream      # Codex：~/.codex/skills/spatial-downstream

# 1. 建环境（只需一次）
conda env create -f "$SKILL_DIR"/envs/spatial-py.yml
Rscript "$SKILL_DIR"/envs/install_r_packages.R      # 只有 P1 / P3 / BANKSY 需要

# 2. 建项目目录
mkdir my_project && cd my_project
cp "$SKILL_DIR"/assets/configs/P1_lowres_scrna.yaml config.yaml
cp "$SKILL_DIR"/assets/samples_template.tsv samples.tsv
cp "$SKILL_DIR"/assets/markers_template.yaml markers.yaml
#    修改 config.yaml（物种、参考路径、半径等）和 samples.tsv（sample、condition、path）

# 3. 运行
bash "$SKILL_DIR"/pipelines/P1_lowres_scrna.sh
```

`samples.tsv` 有三列，用制表符分隔：`sample`、`condition`、`path`。`path` 是 Space Ranger 的 `outs/` 目录、Xenium 输出目录、CosMx 目录，或 CODEX 的 CSV 文件。

每个阶段都会在 `results/` 里写一个 h5ad，所以可以从任意阶段重新开始，手动运行后面的脚本即可。

### 示例：4 个 Visium 样本 + 单细胞参考

假设你有 4 个人肿瘤的 Visium 样本（2 个对照、2 个处理），还有一份同组织的单细胞数据 `ref/scrna.h5ad`，带原始 counts，细胞类型在 `cell_type` 列。你想知道处理之后肿瘤周围的免疫细胞有没有变化。

Visium 是 spot 级数据，又有参考，所以用 P1。可以把这个问题和数据路径直接告诉助手，也可以自己建项目：

```bash
SKILL_DIR=~/.claude/skills/spatial-downstream
mkdir tumor_visium && cd tumor_visium
cp "$SKILL_DIR"/assets/configs/P1_lowres_scrna.yaml config.yaml
cp "$SKILL_DIR"/assets/samples_template.tsv samples.tsv
cp "$SKILL_DIR"/assets/markers_template.yaml markers.yaml
```

`samples.tsv` 每行写一个样本：

```
sample	condition	path
C1	Control	/data/visium/C1/outs
C2	Control	/data/visium/C2/outs
T1	Treated	/data/visium/T1/outs
T2	Treated	/data/visium/T2/outs
```

`config.yaml` 里核对 `species: human`、`scrna_reference: ref/scrna.h5ad`、`scrna_celltype_key: cell_type`，`n_niches` 先用默认的 8，看过 niche 图再调。然后运行 `bash "$SKILL_DIR"/pipelines/P1_lowres_scrna.sh`。跑完先看 `results/qc_violin.png`，过滤阈值不适合你的数据，就改 `preprocess_rna.py` 顶部的参数，从这一步重跑。

用 `results/` 里的这些文件回答问题：

| 问题 | 看哪里 |
|---|---|
| 每个 spot 里有哪些细胞？ | `03_annotated.h5ad` 的 `obsm['cell_props']`；RCTD 报告 `reports/run_rctd.html` |
| 有哪几类微环境（niche），各由什么细胞组成？ | `niche_composition_heatmap.png`、`spatial_niche.png` |
| 每个样本里各 niche 占多少，对照和处理有什么不同？ | `niche_by_sample.csv`、`niche_by_sample.png` |
| 同一个 niche 在不同样本里组成一样吗？ | `niche_composition_by_sample.csv` |
| 哪些细胞类型挨在一起？ | `misty_juxta_mean.png`（各样本的平均） |
| 肿瘤和免疫细胞之间有哪些配受体？ | `lr_local_top.png`、`lr_global_<样本>.csv` |

niche 是 4 个样本一起聚类的，C1 里的 niche 3 和 T2 里的 niche 3 是同一种微环境。RCTD、MISTy 和 Moran's I 是逐样本计算的，表里给的是平均值。每组只有 2 个样本，`niche_by_sample.csv` 只能看趋势，做显著性检验不够可靠，每组 3 个以上样本更可靠。spot 不是单个细胞，所以 P1 不跑邻域富集，细胞共定位看 MISTy 的结果。

### 交互式使用：JupyterLab 笔记本和 R Markdown 报告

`scripts/` 里的每个 Python 脚本，在 `notebooks/` 下都有一个代码相同的 JupyterLab 笔记本，目录结构一样。笔记本开头是说明，接着是设置 `project_dir` 的格子。每个代码格以"步骤 / 输入 / 输出"三行注释开头，后面用一句话说明这一步得到了什么。流程运行的是 `.py`，改了 `.py` 之后运行 `python tools/make_notebooks.py` 重新生成笔记本，两边保持一致。

R 步骤是 R Markdown 报告（`run_rctd.Rmd`、`domain_banksy.Rmd`、`seurat_to_mtx.Rmd`），有编号小节、浮动目录，默认显示代码，每段代码前有说明，有检查用的表格，最后是 `sessionInfo()`。可以在 RStudio 里设好 YAML `params` 中的 `project_dir` 后 Knit，也可以由流程自动渲染到 `results/reports/`。渲染需要 pandoc，RStudio 和 Quarto 都自带；命令行里找不到 pandoc 时，设置 `RSTUDIO_PANDOC`。

## 5. 目录结构

```
spatial-downstream/
├── README.md、README.zh-CN.md   说明文档（英文 / 中文）
├── SKILL.md                     给 Claude Code / Codex 助手的说明
├── VERSION、CHANGELOG.md        版本号和更新记录
├── CITATION.cff、codemeta.json  引用信息和机器可读元数据
├── LICENSE                      MIT 许可证
├── pipelines/                   P1–P7 流程脚本（bash）
├── notebooks/                   JupyterLab 笔记本，每个 Python 脚本一个（自动生成）
├── tools/                       make_notebooks.py：从脚本重新生成笔记本
├── scripts/
│   ├── A_load/                  各平台读入 → 01_raw.h5ad
│   ├── B_preprocess/            质控、标准化、Harmony、Leiden → 02_clustered.h5ad
│   ├── C_annotate/              RCTD、cell2location、标签转移、marker → 03_annotated.h5ad
│   ├── D_domain_niche/          空间域、niche → 04_domain.h5ad、05_niche.h5ad
│   ├── E_neighborhood/          邻域富集、组间差异、MISTy
│   ├── F_communication/         配体-受体分析
│   ├── G_special/               CNV、SVG、热点、NMF、亚型、空间趋势
│   ├── bridge/                  R ↔ Python 数据交换
│   └── Z_provenance/            配置检查、溯源记录
├── docs/                        流程图（SVG）和生成脚本 make_figures.py
├── assets/                      配置预设（configs/P*.yaml）、模板
├── schemas/                     配置的 JSON Schema、AnnData 数据约定
├── envs/                        固定版本的 Python 环境、R 安装脚本、测试过的版本表
├── references/                  选路线指南、方法适用表、多样本说明、单位、R–Python 联合、脚本对照
└── tests/                       模拟数据生成器和端到端测试
```

## 6. R 和 Python 怎么联合

主线对象始终是 h5ad。需要 R 的步骤（RCTD、BANKSY、基于 Seurat 的工具）由 Python 导出普通文件，R 运行后返回以条形码为键的 CSV，Python 再按条形码合并回去：

```
h5ad ──导出（counts.mtx 基因×细胞、genes.txt、barcodes.txt、obs.csv、coords.csv）──▶ R ──result.csv──▶ 按条形码合并 ──▶ h5ad
```

只交换原始整数 counts，两边基因名保持一致，坐标用 µm，合并时一律按条形码，不按行号。完整说明见 [`references/r_python_bridge.md`](references/r_python_bridge.md)。

## 7. 约定

- 单位：所有半径和距离都用 `obsm['spatial_um']`（µm），读入脚本负责把像素、DNB 等平台单位换算成 µm。见 [`references/units_and_params.md`](references/units_and_params.md)。
- 多样本：邻居一律不跨样本。做法是按样本建邻居图（`library_key="sample"`），或先把各样本坐标错开（MISTy、BANKSY），或按样本循环。见 [`references/multi_sample.md`](references/multi_sample.md)。
- 原始数据：`layers['counts']`（RNA）和 `layers['raw']`（蛋白）始终不覆盖，h5ad 保留全部基因。
- spot 数据：用细胞类型比例做 niche 和 MISTy，不用主导类型做邻域富集。
- 靶向 panel：不做 CNV、通路富集、RNA 速度；做配受体之前，先查 panel 覆盖了多少配受体对。
- macOS / Windows：Python 用 spawn 方式开子进程时，进程池会把整个脚本重新执行一遍，所以脚本避开了这类进程池（squidpy 关掉了进度条，CNV 脚本设成 fork）。

## 8. 测试

```bash
bash "$SKILL_DIR"/tests/run_tests.sh /path/to/work_dir
```

`tests/make_test_data.py` 按各平台真实的文件格式生成小型模拟数据（Space Ranger HDF5 + spatial 目录、Visium HD 的 binned parquet、Xenium HDF5 + parquet、CosMx CSV、CODEX CSV），外加一份单细胞参考。测试脚本随后把 P1–P7 从头到尾跑一遍，逐个打印 PASS 或 FAIL。

已测试：7 套流程在模拟数据上全部跑通（macOS arm64，版本见 [`envs/tested_versions.tsv`](envs/tested_versions.tsv)）；Visium、Visium HD、Xenium、CosMx、CODEX 格式的读入；RCTD 的 full 和 doublet 模式；BANKSY 的 R 往返。
尚未测试：真实仪器数据、`load_stereo.py`（stereopy）、`run_cell2location.py`、`commot_flow.py`、`spatial_tendency_soapy.py`、Linux 和 Windows。

## 9. FAIR 说明

| 原则 | 本项目的做法 |
|---|---|
| **F** 可查找 | `README`、`CITATION.cff`、`codemeta.json` 里有描述性名称、关键词和版本号；语义化版本（`VERSION`、`CHANGELOG.md`）。把带标签的版本存档到 Zenodo 可以获得永久标识符（DOI，见第 10 节）。 |
| **A** 可获取 | 脚本和文档都是纯文本；MIT 许可（任何用途都免费，包括商用）；只依赖可从 PyPI、CRAN、Bioconductor 安装的开源软件；不需要登录或商业软件。 |
| **I** 可互操作 | 只用社区标准格式：AnnData `.h5ad`、Matrix Market `.mtx`、CSV/TSV、YAML、JSON Schema。有成文的数据约定（[`schemas/anndata_contract.md`](schemas/anndata_contract.md)）和可校验的配置（[`schemas/config.schema.json`](schemas/config.schema.json)）；基因用 HGNC/MGI symbol，可选 Cell Ontology ID；坐标用国际单位（µm）。 |
| **R** 可重用 | 许可证和引用信息明确；环境版本固定并经过测试；每次运行都会写 `results/provenance/`（配置副本、软件版本、R sessionInfo、输入和输出文件的 SHA-256 校验值）；有基于模拟数据的端到端测试；每个脚本开头写清输入、输出和参数。 |

## 10. 引用和许可证

引用本项目请使用 [`CITATION.cff`](CITATION.cff) 中的信息，同时引用你用到的方法：

| 步骤 | 方法 |
|---|---|
| 预处理 | Scanpy（Wolf et al., *Genome Biol* 2018）；Harmony（Korsunsky et al., *Nat Methods* 2019） |
| 空间统计 | Squidpy（Palla et al., *Nat Methods* 2022） |
| 去卷积 | RCTD（Cable et al., *Nat Biotechnol* 2022）；cell2location（Kleshchevnikov et al., *Nat Biotechnol* 2022） |
| 空间域 | BANKSY（Singhal et al., *Nat Genet* 2024）；CellCharter（Varrone et al., *Nat Genet* 2024） |
| 共定位、配受体 | MISTy（Tanevski et al., *Genome Biol* 2022）；LIANA+（Dimitrov et al., *Nat Cell Biol* 2024）；CellPhoneDB（Efremova et al., *Nat Protoc* 2020）；COMMOT（Cang et al., *Nat Methods* 2023） |
| CNV | inferCNV（Tirosh et al., *Science* 2016）/ infercnvpy |
| 热点 | Getis & Ord（1992），PySAL `esda` |
| 空间趋势 | SOAPy |

获取 DOI：在 GitHub 仓库上打开 Zenodo 集成，发布标签为 `v1.0.0` 的版本，再把 DOI 填回 `CITATION.cff` 和 `codemeta.json`。

许可证：[MIT](LICENSE)。版权所有 © 2026 Yueling Peng。任何人都可以免费使用、复制、修改和再分发，包括商用，只要保留版权和许可声明。
