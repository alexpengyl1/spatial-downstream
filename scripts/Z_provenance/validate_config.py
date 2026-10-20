# 模块 Z · 检查 config.yaml 和样本表（每套流程的第一步）   [Python]
# 作用：在跑几个小时之前，先确认配置字段齐全、取值合法、和所选流程一致，样本路径都存在。
# 依据：schemas/config.schema.json

# %%
# 步骤：导入包
# 输入：无
# 输出：无
from pathlib import Path
import os
import json
import jsonschema
import pandas as pd
import yaml

# %%
# 步骤：读取配置和 schema
# 输入：config.yaml；skill_dir（本 skill 的安装位置）下的 schemas/config.schema.json
# 输出：cfg、schema
skill_dir = Path(os.environ.get("SPATIAL_SKILL_DIR", "~/.claude/skills/spatial-downstream")).expanduser()   # 流程脚本会自动设好；单独运行且装在 Codex 里时改成 ~/.codex/skills/spatial-downstream
cfg = yaml.safe_load(open("config.yaml"))
schema = json.load(open(skill_dir / "schemas" / "config.schema.json"))

# %%
# 步骤：按 schema 检查配置；不合法时直接报错并指出字段
# 输入：cfg、schema
# 输出：无（不合法时抛出 ValidationError）
jsonschema.validate(cfg, schema)

# %%
# 步骤：检查样本表的列名和每个样本路径是否存在
# 输入：samples.tsv
# 输出：samples；打印每个路径是否存在
samples = pd.read_csv(cfg["samples"], sep="\t")
assert list(samples.columns[:3]) == ["sample", "condition", "path"], "samples.tsv 前三列必须是 sample / condition / path"
samples["exists"] = samples["path"].map(lambda p: Path(p).exists())
print(samples)
assert samples["exists"].all(), "有样本路径不存在，见上表"
