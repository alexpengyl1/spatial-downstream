# R 和 Python 联合使用

本 skill 以 Python（scanpy / squidpy）为主线，只有 Python 生态里没有好的替代品时才切到 R：RCTD、BANKSY（可选）、Seurat 对象、inferCNV / CytoTRACE / MultiNicheNet 等 R 包。

## 原则

1. **主线对象只有一个：h5ad。** R 端不保存"另一份完整对象"，只返回结果表（CSV，第一列为条形码）。
2. **交换用最朴素的文件格式。** counts 用 Matrix Market（`.mtx`），基因名、条形码用 txt，元数据和坐标用 csv。两边都能直接读，不依赖版本经常出问题的 SeuratDisk / sceasy。
3. **永远按条形码对齐，不按行号。** R 端过滤（比如 RCTD 丢掉低 UMI 的 spot）之后，行号就对不上了。Python 读回时用 `.reindex(adata.obs_names)` 或 `.loc[...]`。
4. **给 R 的必须是原始整数 counts**（`layers['counts']`），不是 log 值或 scale 后的值。
5. **基因名两边要一致**：都用 symbol，或都用 Ensembl ID；小鼠 / 人不要混用。
6. **坐标统一用 µm**（`obsm['spatial_um']`），R 端的半径、距离参数和 Python 端含义相同。
7. **矩阵方向**：AnnData 是 细胞 x 基因，R / Seurat 是 基因 x 细胞。导出时转置（`.T`），读回时再转置。

## 方式一（默认）：文件交换

```
Python                                   R                                       Python
03_annotated.h5ad ──h5ad_to_mtx.py──▶ exchange/ ──domain_banksy.Rmd──▶ banksy_domain.csv ──import_obs_csv.py──▶ 04_domain.h5ad
02_clustered.h5ad ──export_for_rctd.py──▶ rctd/ ──run_rctd.Rmd──▶ weights.csv ──import_rctd.py──▶ 03_annotated.h5ad
Seurat .rds ──seurat_to_mtx.Rmd──▶ exchange_seurat/ ──mtx_to_h5ad.py──▶ scrna.h5ad（作为单细胞参考）
```

R 端读入的固定写法：

```r
counts <- as(Matrix::readMM("counts.mtx"), "CsparseMatrix")   # 基因 x 细胞
rownames(counts) <- readLines("genes.txt")
colnames(counts) <- readLines("barcodes.txt")
obs <- read.csv("obs.csv", row.names = 1)
coords <- read.csv("coords.csv", row.names = 1)                 # x, y（µm）
```

R 端写回的固定写法（一行一个细胞，行名为条形码）：

```r
write.csv(data.frame(domain = clusters, row.names = colnames(counts)), "result.csv")
```

想把 R 的结果合并回 h5ad 时，改 `bridge/import_obs_csv.py` 顶部的三个路径就行。

## 方式二：R 里直接读写 h5ad（zellkonverter）

适合在 R 里交互式分析，不想经过中间文件：

```r
# 步骤：把 h5ad 读成 SingleCellExperiment（X 进 assay "X"，layers 各成一个 assay）
# 输入：results/03_annotated.h5ad
# 输出：sce
sce <- zellkonverter::readH5AD("results/03_annotated.h5ad")
SummarizedExperiment::assayNames(sce)          # 应该能看到 "X" 和 "counts"

# 步骤：转成 Seurat（如果后续要用 Seurat 函数）
# 输入：sce
# 输出：seu
seu <- Seurat::as.Seurat(sce, counts = "counts", data = "X")

# 步骤：结果写回 h5ad（R 端新加的 colData 列会进 obs）
# 输入：sce
# 输出：results/03_from_r.h5ad
zellkonverter::writeH5AD(sce, "results/03_from_r.h5ad")
```

注意：zellkonverter 第一次运行会用 basilisk 自动装一个 Python 环境，比较慢；`obsm` 里的 DataFrame（比如 `cell_props`）会变成矩阵，列名可能丢失。所以流程化运行时仍然用方式一。

## 方式三：Seurat 空间对象 → Python

如果前面的分析是在 Seurat 里做的（课程第4、13课的 RCTD 就是读 Seurat rds），用 `bridge/seurat_to_mtx.Rmd`：

- `assay` 改成空间对象的 assay 名（Visium 为 `Spatial`，Visium HD 为 `Spatial.016um` 或 `Spatial.008um`）
- 取消注释写坐标的那一行；`GetTissueCoordinates()` 给的是全分辨率像素坐标，在 Python 里乘以 µm/像素（Visium：55 / spot_diameter_fullres；HD：scalefactors 里的 microns_per_pixel）
- 再用 `bridge/mtx_to_h5ad.py` 组装，并在其中加上 `adata.obsm["spatial_um"] = ...`

## 不推荐的方式

- **reticulate / rpy2 在同一个进程里来回调用**：环境冲突难排查，出了错不知道是哪边的问题。
- **SeuratDisk 的 h5Seurat ↔ h5ad**：对 Seurat v5 的多 layer 结构支持不好，常常丢失 counts。
- **用 CSV 存完整表达矩阵**：几万个细胞时文件就有几个 GB，读写都慢。

## 何时必须用 R

| 需求 | R 包 | 原因 | 课程 |
|---|---|---|---|
| spot 去卷积（无 GPU） | spacexr（RCTD） | 速度快，doublet 模式适合 HD | 4、5、13 |
| 空间域（官方实现） | Banksy | Bioconductor 版维护最好 | 12 |
| 多样本多分组的细胞通讯 | multinichenetr | 只有 R 版 | 2 |
| 细胞干性 | CytoTRACE | 原版只有 R | 2 |
| 空间模式基因 + 配受体共定位 | SpaGene | 只有 R 版（Python 用 svg_moran.py + lr_spatial_liana.py 代替） | 7 |

这些都按"Python 导出 → R 运行 → CSV 回传 → Python 合并"的流程处理。
