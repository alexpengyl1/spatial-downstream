# 模块 Z · 记录溯源信息（每套流程的最后一步）   [Python]
# 作用：把"用什么数据、什么参数、什么软件版本、产出了什么"写进 results/provenance/，让结果可以被复查和复现（FAIR 中的 R）。
# 运行：在项目目录下 python record_provenance.py（流程编号取自 config.yaml 的 pipeline 字段）

# %%
# 步骤：导入包
# 输入：无
# 输出：无
from datetime import datetime
from importlib.metadata import distributions
from pathlib import Path
import os
import hashlib
import platform
import shutil
import subprocess
import pandas as pd
import yaml

# %%
# 步骤：读取配置，建立溯源目录
# 输入：config.yaml；skill_dir（本 skill 的安装位置）
# 输出：cfg、prov_dir、pipeline_id、skill_dir
cfg = yaml.safe_load(open("config.yaml"))
prov_dir = Path(cfg["out_dir"]) / "provenance"
prov_dir.mkdir(parents=True, exist_ok=True)
pipeline_id = cfg["pipeline"]
skill_dir = Path(os.environ.get("SPATIAL_SKILL_DIR", "~/.claude/skills/spatial-downstream")).expanduser()   # 流程脚本会自动设好；单独运行且装在 Codex 里时改成 ~/.codex/skills/spatial-downstream
max_hash_mb = 100          # 大于 100MB 的文件只记大小，不算 SHA-256（太慢）

# %%
# 步骤：复制本次运行用到的配置文件
# 输入：config.yaml、samples.tsv、markers.yaml（如果有）
# 输出：provenance/ 下的同名文件
for name in ["config.yaml", cfg["samples"], cfg.get("markers", "")]:
    if name and Path(name).exists():
        shutil.copy(name, prov_dir / Path(name).name)

# %%
# 步骤：记录运行信息
# 输入：无
# 输出：provenance/run_info.yaml
run_info = {
    "pipeline": pipeline_id,
    "skill_version": (skill_dir / "VERSION").read_text().strip(),
    "finished_at": datetime.now().isoformat(timespec="seconds"),
    "os": platform.platform(),
    "python": platform.python_version(),
}
yaml.safe_dump(run_info, open(prov_dir / "run_info.yaml", "w"), allow_unicode=True)
run_info

# %%
# 步骤：记录 Python 包版本
# 输入：envs/tested_versions.tsv 里列出的 Python 包
# 输出：provenance/python_versions.tsv（包名、本次版本、测试过的版本）
installed = {d.metadata["Name"].lower(): d.version for d in distributions()}
tested = pd.read_csv(skill_dir / "envs" / "tested_versions.tsv", sep="\t")
tested = tested[(tested["language"] == "Python") & (tested["package"] != "python")]
tested["installed"] = tested["package"].str.lower().map(installed).fillna("not installed")
tested[["package", "installed", "version"]].rename(columns={"version": "tested"}).to_csv(prov_dir / "python_versions.tsv", sep="\t", index=False)
tested[["package", "installed", "version"]]

# %%
# 步骤：如果本机有 R，记录 R 的 sessionInfo
# 输入：Rscript
# 输出：provenance/r_session_info.txt
if shutil.which("Rscript"):
    r_code = 'for (p in c("spacexr", "Banksy")) if (requireNamespace(p, quietly = TRUE)) cat(p, as.character(packageVersion(p)), "\\n"); print(sessionInfo())'
    result = subprocess.run(["Rscript", "-e", r_code], capture_output=True, text=True)
    (prov_dir / "r_session_info.txt").write_text(result.stdout + result.stderr)

# %%
# 步骤：列出输入文件（各样本路径下的所有文件）和输出文件（out_dir 下除 provenance/ 以外的所有文件）
# 输入：samples.tsv、out_dir
# 输出：file_lists（清单文件名 -> 文件列表）
samples = pd.read_csv(cfg["samples"], sep="\t")
file_lists = {
    "input_files.tsv": [f for p in samples["path"] for f in ([Path(p)] if Path(p).is_file() else Path(p).rglob("*")) if f.is_file()],
    "output_files.tsv": [f for f in Path(cfg["out_dir"]).rglob("*") if f.is_file() and prov_dir not in f.parents],
}
{name: len(files) for name, files in file_lists.items()}

# %%
# 步骤：每个文件记录路径、大小、修改时间和 SHA-256（大于 max_hash_mb 的文件不算校验值）
# 输入：file_lists
# 输出：provenance/input_files.tsv、provenance/output_files.tsv
for name, files in file_lists.items():
    table = pd.DataFrame({"path": [str(f) for f in files],
                          "size_mb": [round(f.stat().st_size / 1e6, 3) for f in files],
                          "modified": [datetime.fromtimestamp(f.stat().st_mtime).isoformat(timespec="seconds") for f in files]})
    table["sha256"] = [hashlib.sha256(f.read_bytes()).hexdigest() if size <= max_hash_mb else "" for f, size in zip(files, table["size_mb"])]
    table.to_csv(prov_dir / name, sep="\t", index=False)
