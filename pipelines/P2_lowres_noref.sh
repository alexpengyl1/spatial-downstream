#!/usr/bin/env bash
# P2 · 低精度空转，无单细胞参考 / Low-resolution spatial, no scRNA-seq reference
# 平台 / platforms：10X Visium、Stereo-seq bin50 及以上
# 注释 / annotation：marker 打分（obsm['cell_scores']，不是比例，只能近似细胞组成）；以表达聚类和空间域为主线
# 配置 / config：assets/configs/P2_lowres_noref.yaml（composition_key: cell_scores）+ markers.yaml
# 用法 / usage：在项目目录下运行 bash P2_lowres_noref.sh
set -e

SKILL_DIR="$(cd "$(dirname "$0")/.." && pwd)"   # 本 skill 的安装位置（~/.claude/skills/ 或 ~/.codex/skills/ 都行），由本脚本所在位置推出
export SPATIAL_SKILL_DIR="$SKILL_DIR"         # 告诉 Python 脚本 schema 等文件在哪
SKILL="$SKILL_DIR/scripts"
# PY 和 R_CMD 是带空格的命令，下面调用时故意不加引号，让 shell 拆成多个词
PY="${PY:-conda run --no-capture-output -n spatial-py python}"   # 主 Python 环境 / main Python env（可用环境变量覆盖）
LOADER="${LOADER:-load_visium.py}"      # Stereo-seq 改成 load_stereo.py（stereo 环境）

# Z：检查配置   [Python]
$PY "$SKILL/Z_provenance/validate_config.py"

# A + B：读入、质控、聚类   [Python]
$PY "$SKILL/A_load/$LOADER"
$PY "$SKILL/B_preprocess/preprocess_rna.py"

# C：marker 打分   [Python]
$PY "$SKILL/C_annotate/marker_annotation.py"

# D：空间域 + niche（基于 marker 分数）   [Python]
$PY "$SKILL/D_domain_niche/domain_cellcharter.py"
$PY "$SKILL/D_domain_niche/niche_composition.py"

# E：共定位（MISTy，基于 marker 分数）   [Python]
$PY "$SKILL/E_neighborhood/misty_colocalization.py"

# F：空间配受体   [Python]
$PY "$SKILL/F_communication/lr_spatial_liana.py"
# 可选 / optional: conda run -n commot python "$SKILL/F_communication/commot_flow.py"

# G：个性化分析   [Python]
$PY "$SKILL/G_special/svg_moran.py"
$PY "$SKILL/G_special/nmf_programs.py"
# $PY "$SKILL/G_special/cnv_infercnvpy.py"      # 肿瘤：用免疫细胞得分高的类型作参考
# $PY "$SKILL/G_special/hotspot_getis.py"

# Z：记录溯源信息   [Python]
$PY "$SKILL/Z_provenance/record_provenance.py"
