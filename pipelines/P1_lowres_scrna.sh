#!/usr/bin/env bash
# P1 · 低精度空转 + 单细胞参考 / Low-resolution spatial + scRNA-seq reference
# 平台 / platforms：10X Visium（v1/v2/CytAssist）、Stereo-seq bin50 及以上
# 注释 / annotation：RCTD full 模式（Python → R → Python）；有 GPU 时可换 cell2location
# 配置 / config：assets/configs/P1_lowres_scrna.yaml，复制到项目目录并改名为 config.yaml
# 用法 / usage：在项目目录下运行 bash P1_lowres_scrna.sh
set -e

SKILL_DIR="$(cd "$(dirname "$0")/.." && pwd)"   # 本 skill 的安装位置（~/.claude/skills/ 或 ~/.codex/skills/ 都行），由本脚本所在位置推出
export SPATIAL_SKILL_DIR="$SKILL_DIR"         # 告诉 Python 脚本 schema 等文件在哪
SKILL="$SKILL_DIR/scripts"
# PY 和 R_CMD 是带空格的命令，下面调用时故意不加引号，让 shell 拆成多个词
PY="${PY:-conda run --no-capture-output -n spatial-py python}"   # 主 Python 环境 / main Python env（可用环境变量覆盖）
R_CMD="${R_CMD:-Rscript}"   # 装了 spacexr 的 R；装在 conda 环境里时改成 "conda run -n 环境名 Rscript"
# R 步骤是 R Markdown 报告，渲染需要 pandoc：RStudio / Quarto 自带；命令行找不到时用 RSTUDIO_PANDOC 指向它（下面是 macOS + RStudio 的位置）
export LANG="${LANG:-en_US.UTF-8}"     # 报告里有中文，R 需要 UTF-8 区域设置
export RSTUDIO_PANDOC="${RSTUDIO_PANDOC:-/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64}"
LOADER="${LOADER:-load_visium.py}"      # Stereo-seq 改成 load_stereo.py（需要 stereo 环境，可先单独运行这一步）

# Z：检查配置   [Python]
$PY "$SKILL/Z_provenance/validate_config.py"

# A + B：读入、质控、聚类   [Python]
$PY "$SKILL/A_load/$LOADER"
$PY "$SKILL/B_preprocess/preprocess_rna.py"

# C：去卷积   [Python → R → Python]
$PY "$SKILL/C_annotate/prepare_reference.py"
$PY "$SKILL/C_annotate/export_for_rctd.py"
# 先把项目目录存进变量：render() 的参数是惰性求值的，到它切换到 Rmd 所在目录之后才会执行 getwd()
$R_CMD -e "project <- getwd(); rmarkdown::render('$SKILL/C_annotate/run_rctd.Rmd', params = list(project_dir = project), output_dir = file.path(project, yaml::read_yaml('config.yaml')\$out_dir, 'reports'), intermediates_dir = tempdir())"
$PY "$SKILL/C_annotate/import_rctd.py"
# 替代方案（cell2location 环境，GPU）/ alternative: conda run -n cell2location python "$SKILL/C_annotate/run_cell2location.py"

# D：空间域 + 细胞 niche   [Python]
$PY "$SKILL/D_domain_niche/domain_cellcharter.py"
$PY "$SKILL/D_domain_niche/niche_composition.py"

# E：细胞类型共定位（spot 数据用 MISTy）   [Python]
$PY "$SKILL/E_neighborhood/misty_colocalization.py"

# F：空间配受体   [Python]
$PY "$SKILL/F_communication/lr_spatial_liana.py"
# 可选，信号流向（commot 环境，较慢）/ optional: conda run -n commot python "$SKILL/F_communication/commot_flow.py"

# G：个性化分析   [Python]
$PY "$SKILL/G_special/svg_moran.py"
$PY "$SKILL/G_special/nmf_programs.py"
# 可选（先改脚本顶部参数或 config）/ optional, edit parameters first:
# $PY "$SKILL/G_special/cnv_infercnvpy.py"      # 肿瘤：config 里填 gtf、cnv_reference
# $PY "$SKILL/G_special/hotspot_getis.py"       # 改 feature_a / feature_b

# Z：记录溯源信息   [Python]
$PY "$SKILL/Z_provenance/record_provenance.py"
