#!/usr/bin/env bash
# P7 · 空间蛋白组 / Spatial proteomics
# 平台 / platforms：CODEX（PhenoCycler）、IMC、mIHC、MIBI
# 注释 / annotation：marker 打分 + 人工核对 marker_matrixplot.png（单细胞 RNA 参考一般不适用）
# 不做 / not applicable：去卷积、CNV、配受体通讯（抗体 panel 里基本没有配受体）
# 配置 / config：assets/configs/P7_protein.yaml（pixel_size_um 按显微镜填写）+ markers.yaml（写抗体名）
# 用法 / usage：在项目目录下运行 bash P7_protein.sh
set -e

SKILL_DIR="$(cd "$(dirname "$0")/.." && pwd)"   # 本 skill 的安装位置（~/.claude/skills/ 或 ~/.codex/skills/ 都行），由本脚本所在位置推出
export SPATIAL_SKILL_DIR="$SKILL_DIR"         # 告诉 Python 脚本 schema 等文件在哪
SKILL="$SKILL_DIR/scripts"
# PY 和 R_CMD 是带空格的命令，下面调用时故意不加引号，让 shell 拆成多个词
PY="${PY:-conda run --no-capture-output -n spatial-py python}"   # 主 Python 环境 / main Python env（可用环境变量覆盖）
HAS_GROUPS="${HAS_GROUPS:-yes}"

# Z：检查配置   [Python]
$PY "$SKILL/Z_provenance/validate_config.py"

# A + B：读入、arcsinh 转换、聚类   [Python]
$PY "$SKILL/A_load/load_codex.py"
$PY "$SKILL/B_preprocess/preprocess_protein.py"

# C：marker 注释   [Python]
$PY "$SKILL/C_annotate/marker_annotation.py"

# D：空间域 + CN（cellular neighborhood）   [Python]
$PY "$SKILL/D_domain_niche/domain_cellcharter.py"
$PY "$SKILL/D_domain_niche/niche_composition.py"

# E：邻域富集、组间差异、MISTy   [Python]
$PY "$SKILL/E_neighborhood/nhood_enrichment.py"
if [ "$HAS_GROUPS" = "yes" ]; then
    $PY "$SKILL/E_neighborhood/diff_nhood.py"
fi
$PY "$SKILL/E_neighborhood/misty_colocalization.py"

# G：个性化分析   [Python]
$PY "$SKILL/G_special/svg_moran.py"
# $PY "$SKILL/G_special/nmf_subtype.py"        # 改 target_cell_type

# Z：记录溯源信息   [Python]
$PY "$SKILL/Z_provenance/record_provenance.py"
