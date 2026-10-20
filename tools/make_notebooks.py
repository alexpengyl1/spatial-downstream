# 把 scripts/ 下的每个 .py 生成一个 JupyterLab 笔记本（notebooks/ 下同名 .ipynb）
# 运行：python tools/make_notebooks.py
# 说明：流程（pipelines/*.sh）运行的是 .py；笔记本是同一份代码，加上标题、每步的小节和结果说明，方便在 JupyterLab 里逐格运行。
#      改了 .py 之后重新运行本脚本，笔记本就会同步，不需要手动维护两份代码。
#      笔记本的写法遵循 readable-coding：每个代码格开头是"步骤 / 输入 / 输出"三行注释，代码格之后用一句话说明得到了什么。

# %%
# 步骤：导入包，确定输入和输出目录
# 输入：无
# 输出：skill_dir、scripts_dir、notebooks_dir
from pathlib import Path
import json
import re

skill_dir = Path(__file__).resolve().parents[1]
scripts_dir = skill_dir / "scripts"
notebooks_dir = skill_dir / "notebooks"

# %%
# 步骤：笔记本第一个代码格：切换到项目目录（脚本默认在项目目录下运行，笔记本则在自己所在的目录启动）
# 输入：无
# 输出：project_cell（代码格的源码）
project_cell = '''# 步骤：切换到项目目录（里面有 config.yaml、samples.tsv 和 results/）
# 输入：project_dir（改成自己的项目目录）
# 输出：无（之后的代码都以项目目录为当前目录）
import os
project_dir = "."      # 例如 "/data/my_visium_project"
os.chdir(project_dir)
os.getcwd()'''

# %%
# 步骤：生成 Jupyter 格子（nbformat 4 的 JSON 结构）
# 输入：格子类型和源码
# 输出：markdown()、code() 返回一个格子
def markdown(text):
    return {"cell_type": "markdown", "metadata": {}, "source": text}

def code(text):
    return {"cell_type": "code", "execution_count": None, "metadata": {}, "outputs": [], "source": text}

# %%
# 步骤：逐个脚本转换：开头的注释变成标题和说明；每个 "# %%" 代码块前加小节标题，后加一句结果说明
# 输入：scripts/**/*.py
# 输出：notebooks/<模块>/<脚本名>.ipynb
for script in sorted(scripts_dir.glob("*/*.py")):
    text = script.read_text()
    intro, *blocks = re.split(r"^# %%\n", text, flags=re.M)
    intro_lines = []
    for line in intro.strip().splitlines():
        if line.startswith("#  ") and intro_lines:          # 缩进的注释行是上一行的续行
            intro_lines[-1] += " " + line.lstrip("# ").rstrip()
        elif line.strip("# "):
            intro_lines.append(line.lstrip("# ").rstrip())
    cells = [markdown(f"# {intro_lines[0]}\n\n" + "\n".join(f"- {line}" for line in intro_lines[1:]) +
                      f"\n\n**对应的流程脚本**：`scripts/{script.parent.name}/{script.name}`（代码相同，流程运行的是 .py）\n\n"
                      "**用法**：先把下面第一个代码格里的 `project_dir` 改成项目目录，然后从上到下逐格运行。"),
             code(project_cell)]
    for i, block in enumerate(blocks, start=1):
        block = block.strip("\n")
        step = re.search(r"^# 步骤：(.+)$", block, flags=re.M).group(1)
        output = re.search(r"^# 输出：(.+)$", block, flags=re.M).group(1)
        cells += [markdown(f"## {i}. {step}"), code(block)]
        if output.strip() != "无":
            cells.append(markdown(f"这一步得到：{output}"))
    notebook = {"cells": cells, "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
                                              "language_info": {"name": "python"}}, "nbformat": 4, "nbformat_minor": 5}
    out = notebooks_dir / script.parent.name / (script.stem + ".ipynb")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(notebook, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")

# %%
# 步骤：列出生成的笔记本
# 输入：notebooks/
# 输出：打印的文件清单
print(len(list(notebooks_dir.glob("*/*.ipynb"))), "notebooks")
