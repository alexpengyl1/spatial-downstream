#!/usr/bin/env bash
# P3 · 高精度测序型空转 + 单细胞参考 / High-resolution sequencing-based spatial + scRNA-seq reference
# 平台 / platforms：Visium HD（16µm / 8µm bin）、Stereo-seq bin20 及以下 / cellbin、Slide-seq
# 注释 / annotation：RCTD doublet 模式（Python → R → Python）；用细胞分割结果时可改用 label_transfer.py
# 配置 / config：assets/configs/P3_hd_scrna.yaml（rctd_mode: doublet）
# 用法 / usage：在项目目录下运行 bash P3_hd_scrna.sh
# 耗时 / runtime：RCTD doublet 1 万个 bin、4 核约 3.5 分钟；几十万个 bin 需数小时
set -e

SKILL_DIR="$(cd "$(dirname "$0")/.." && pwd)"   # 本 skill 的安装位置（~/.claude/skills/ 或 ~/.codex/skills/ 都行），由本脚本所在位置推出
export SPATIAL_SKILL_DIR="$SKILL_DIR"         # 告诉 Python 脚本 schema 等文件在哪
SKILL="$SKILL_DIR/scripts"
# PY 和 R_CMD 是带空格的命令，下面调用时故意不加引号，让 shell 拆成多个词
PY="${PY:-conda run --no-capture-output -n spatial-py python}"   # 主 Python 环境 / main Python env（可用环境变量覆盖）
R_CMD="${R_CMD:-Rscript}"
# R 步骤是 R Markdown 报告，渲染需要 pandoc：RStudio / Quarto 自带；命令行找不到时用 RSTUDIO_PANDOC 指向它（下面是 macOS + RStudio 的位置）
export LANG="${LANG:-en_US.UTF-8}"     # 报告里有中文，R 需要 UTF-8 区域设置
export RSTUDIO_PANDOC="${RSTUDIO_PANDOC:-/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64}"
LOADER="${LOADER:-load_visium_hd.py}"   # Stereo-seq 改成 load_stereo.py（stereo 环境）
HAS_GROUPS="${HAS_GROUPS:-yes}"             # samples.tsv 里有 ≥ 2 个分组、每组 ≥ 2 个样本时写 yes

# Z：检查配置   [Python]
$PY "$SKILL/Z_provenance/validate_config.py"

# A + B：读入、质控、聚类   [Python]
$PY "$SKILL/A_load/$LOADER"
$PY "$SKILL/B_preprocess/preprocess_rna.py"

# C：RCTD doublet 注释   [Python → R → Python]
$PY "$SKILL/C_annotate/prepare_reference.py"
$PY "$SKILL/C_annotate/export_for_rctd.py"
# 先把项目目录存进变量：render() 的参数是惰性求值的，到它切换到 Rmd 所在目录之后才会执行 getwd()
$R_CMD -e "project <- getwd(); rmarkdown::render('$SKILL/C_annotate/run_rctd.Rmd', params = list(project_dir = project), output_dir = file.path(project, yaml::read_yaml('config.yaml')\$out_dir, 'reports'), intermediates_dir = tempdir())"
$PY "$SKILL/C_annotate/import_rctd.py"

# D：空间域 + 细胞 niche   [Python]
$PY "$SKILL/D_domain_niche/domain_cellcharter.py"
$PY "$SKILL/D_domain_niche/niche_composition.py"

# E：邻域富集 + 组间差异   [Python]
$PY "$SKILL/E_neighborhood/nhood_enrichment.py"
if [ "$HAS_GROUPS" = "yes" ]; then
    $PY "$SKILL/E_neighborhood/diff_nhood.py"
fi

# F：空间配受体   [Python]
$PY "$SKILL/F_communication/lr_spatial_liana.py"
# 先在脚本顶部改 target_cell_type / edit target_cell_type first:
# $PY "$SKILL/F_communication/ligrec_neighborhood.py"

# G：个性化分析   [Python]
$PY "$SKILL/G_special/svg_moran.py"
$PY "$SKILL/G_special/nmf_programs.py"
# $PY "$SKILL/G_special/nmf_subtype.py"        # 改 target_cell_type
# $PY "$SKILL/G_special/cnv_infercnvpy.py"     # 肿瘤：填 gtf、cnv_reference

# Z：记录溯源信息   [Python]
$PY "$SKILL/Z_provenance/record_provenance.py"
