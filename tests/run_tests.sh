#!/usr/bin/env bash
# 用模拟数据把 7 套流程完整跑一遍 / Run all 7 pipelines end-to-end on simulated data
# 用法 / usage：bash run_tests.sh [工作目录]   （默认 ./spatial_downstream_tests）
# 可用环境变量覆盖 / overridable: PY（Python 命令）、R_CMD（R 命令）
# 结果 / result：每套流程打印 PASS 或 FAIL，日志在 <工作目录>/logs/
set -u

TESTS_DIR="$(cd "$(dirname "$0")" && pwd)"
PIPELINES_DIR="$TESTS_DIR/../pipelines"
WORK_DIR="${1:-./spatial_downstream_tests}"
export PY="${PY:-conda run --no-capture-output -n spatial-py python}"
export R_CMD="${R_CMD:-Rscript}"
export MPLBACKEND=Agg
export LANG="${LANG:-en_US.UTF-8}"
export RSTUDIO_PANDOC="${RSTUDIO_PANDOC:-/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64}"   # R Markdown 渲染用的 pandoc

# 步骤：生成模拟数据和 7 个测试项目
mkdir -p "$WORK_DIR/logs"
$PY "$TESTS_DIR/make_test_data.py" "$WORK_DIR"

# 步骤：逐个流程运行（P6 用 CosMx 读入脚本，顺便测试它）
for pipeline in P1_lowres_scrna P2_lowres_noref P3_hd_scrna P4_hd_noref P5_imaging_scrna P6_imaging_noref P7_protein; do
    loader=""
    if [ "$pipeline" = "P6_imaging_noref" ]; then loader="load_cosmx.py"; fi
    start=$(date +%s)
    if (cd "$WORK_DIR/projects/$pipeline" && LOADER="${loader:-}" bash "$PIPELINES_DIR/$pipeline.sh") > "$WORK_DIR/logs/$pipeline.log" 2>&1; then
        echo "PASS  $pipeline  ($(( $(date +%s) - start ))s)"
    else
        echo "FAIL  $pipeline  （见 / see $WORK_DIR/logs/$pipeline.log）"
    fi
done
