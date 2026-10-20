#!/usr/bin/env bash
# P5 · 成像型单细胞空转 + 单细胞参考 / Imaging-based single-cell spatial + scRNA-seq reference
# 平台 / platforms：10X Xenium、NanoString CosMx、Vizgen MERSCOPE（靶向基因 panel）
# 注释 / annotation：Harmony + KNN 标签转移（看 cell_type_score 判断置信度）
# 不做 / not applicable：CNV、通路富集、COMMOT（panel 基因太少）
# 配置 / config：assets/configs/P5_imaging_scrna.yaml
# 用法 / usage：在项目目录下运行 bash P5_imaging_scrna.sh
set -e

SKILL_DIR="$(cd "$(dirname "$0")/.." && pwd)"   # 本 skill 的安装位置（~/.claude/skills/ 或 ~/.codex/skills/ 都行），由本脚本所在位置推出
export SPATIAL_SKILL_DIR="$SKILL_DIR"         # 告诉 Python 脚本 schema 等文件在哪
SKILL="$SKILL_DIR/scripts"
# PY 和 R_CMD 是带空格的命令，下面调用时故意不加引号，让 shell 拆成多个词
PY="${PY:-conda run --no-capture-output -n spatial-py python}"   # 主 Python 环境 / main Python env（可用环境变量覆盖）
LOADER="${LOADER:-load_xenium.py}"      # CosMx 改成 load_cosmx.py
HAS_GROUPS="${HAS_GROUPS:-yes}"

# Z：检查配置   [Python]
$PY "$SKILL/Z_provenance/validate_config.py"

# A + B：读入、质控、聚类   [Python]
$PY "$SKILL/A_load/$LOADER"
$PY "$SKILL/B_preprocess/preprocess_rna.py"

# C：标签转移   [Python]
$PY "$SKILL/C_annotate/prepare_reference.py"
$PY "$SKILL/C_annotate/label_transfer.py"

# D：空间域 + RCN   [Python]
$PY "$SKILL/D_domain_niche/domain_cellcharter.py"
$PY "$SKILL/D_domain_niche/niche_composition.py"

# E：邻域富集 + 组间差异   [Python]
$PY "$SKILL/E_neighborhood/nhood_enrichment.py"
if [ "$HAS_GROUPS" = "yes" ]; then
    $PY "$SKILL/E_neighborhood/diff_nhood.py"
fi

# F：空间配受体（先看 lr_global_*.csv 里 panel 覆盖了多少配受体对）   [Python]
$PY "$SKILL/F_communication/lr_spatial_liana.py"
# $PY "$SKILL/F_communication/ligrec_neighborhood.py"   # 改 target_cell_type

# G：个性化分析   [Python]
$PY "$SKILL/G_special/svg_moran.py"
# $PY "$SKILL/G_special/nmf_subtype.py"        # 改 target_cell_type

# Z：记录溯源信息   [Python]
$PY "$SKILL/Z_provenance/record_provenance.py"
