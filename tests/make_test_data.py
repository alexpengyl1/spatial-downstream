# 测试数据生成器：按各平台真实的输出格式生成小型模拟数据，并建好 7 个测试项目
# Test-data generator: writes small simulated datasets in each platform's real output format and sets up 7 test projects
# 运行 / usage：python make_test_data.py <工作目录 / work dir>
# 模拟设计 / design：5 种细胞（Tumor 在中心，外面一圈 Fibroblasts，最外围 Macrophage，T cells 和 Endothelial 随机散布）；
#                   Treated 组的肿瘤核心更大。每种细胞有 12 个真实 marker 基因（多为配受体），再加 240 个背景基因。

# %%
# 步骤：导入包，设定工作目录
# 输入：命令行参数（工作目录）
# 输出：work_dir、skill_dir、rng
from pathlib import Path
import json
import sys
import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
import yaml

work_dir = Path(sys.argv[1]).resolve()
skill_dir = Path(__file__).resolve().parents[1]
rng = np.random.default_rng(0)

# %%
# 步骤：定义细胞类型、marker 基因和表达特征
# 输入：无
# 输出：cell_types、genes、signature（细胞类型 x 基因的期望表达）
markers = {
    "Tumor": ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "MUC1", "EGFR", "ERBB2", "CD44", "VEGFA", "MKI67", "TOP2A"],
    "Fibroblasts": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA", "PDGFRB", "FAP", "CXCL12", "TGFBR2", "FN1", "POSTN"],
    "Macrophage": ["CD68", "CD163", "CD14", "LYZ", "C1QA", "C1QB", "SPP1", "TGFB1", "CSF1R", "MRC1", "APOE", "CD74"],
    "T cells": ["CD3D", "CD3E", "CD2", "CD8A", "CD4", "IL7R", "CCL5", "CXCR4", "GZMB", "PTPRC", "CD69", "IFNG"],
    "Endothelial": ["PECAM1", "VWF", "CDH5", "KDR", "FLT1", "ENG", "EMCN", "PLVAP", "ACKR1", "SELE", "ICAM1", "NOTCH1"],
}
cell_types = list(markers)
genes = [g for gs in markers.values() for g in gs] + [f"GENE{i}" for i in range(240)] + ["MT-CO1", "MT-ND1"]
signature = np.full((len(cell_types), len(genes)), 0.3)
for i, cell_type in enumerate(cell_types):
    signature[i, [genes.index(g) for g in markers[cell_type]]] = 4.0

# %%
# 步骤：定义空间布局和计数生成函数
# 输入：坐标（µm）、组别
# 输出：layout() 返回每个点的细胞类型编号；simulate_counts() 返回计数矩阵
def layout(xy, treated, size):
    r = np.hypot(xy[:, 0] - size / 2, xy[:, 1] - size / 2) / size
    core = 0.22 if treated else 0.15
    label = np.where(r < core, 0, np.where(r < core + 0.12, 1, 2))
    u = rng.random(len(xy))
    label = np.where(u < 0.15, 3, label)
    return np.where((u >= 0.15) & (u < 0.25), 4, label)

def simulate_counts(labels_per_point, scale):
    expected = np.stack([signature[labels].sum(axis=0) for labels in labels_per_point]) * scale
    return rng.poisson(expected).astype(np.int32)

# %%
# 步骤：10X HDF5 写入函数（Space Ranger / Xenium 都用这个格式）
# 输入：counts（点 x 基因）、条形码、基因名
# 输出：.h5 文件
def write_10x_h5(path, counts, barcodes):
    matrix = sparse.csc_matrix(counts.T)
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as f:
        group = f.create_group("matrix")
        group["barcodes"] = np.array(barcodes, dtype="S")
        group["data"], group["indices"], group["indptr"] = matrix.data, matrix.indices, matrix.indptr
        group["shape"] = np.array(matrix.shape, dtype="int32")
        features = group.create_group("features")
        features["id"] = np.array([f"ENSG{i:011d}" for i in range(len(genes))], dtype="S")
        features["name"] = np.array(genes, dtype="S")
        features["feature_type"] = np.array(["Gene Expression"] * len(genes), dtype="S")
        features["genome"] = np.array(["GRCh38"] * len(genes), dtype="S")
        features["_all_tag_keys"] = np.array([b"genome"])

# %%
# 步骤：写 Visium 样本（55µm spot，100µm 间距，六边形排列；每个 spot 混合 8 个细胞）
# 输入：样本名、组别
# 输出：<样本>/outs/ 下的 filtered_feature_bc_matrix.h5 和 spatial/
def write_visium(out, treated):
    rows, cols = np.meshgrid(np.arange(12), np.arange(12), indexing="ij")
    xy = np.c_[(cols * 100 + (rows % 2) * 50).ravel(), (rows * 86.6).ravel()].astype(float)
    labels = [layout(xy + rng.normal(0, 30, xy.shape), treated, 1150) for _ in range(8)]
    counts = simulate_counts(np.stack(labels, axis=1), 2.0)      # 约 2000 UMI / spot，接近真实 Visium
    barcodes = [f"SPOT{i:04d}-1" for i in range(len(xy))]
    write_10x_h5(out / "filtered_feature_bc_matrix.h5", counts, barcodes)
    um_per_px = 0.5
    positions = pd.DataFrame({"barcode": barcodes, "in_tissue": 1, "array_row": rows.ravel(), "array_col": cols.ravel(),
                              "pxl_row_in_fullres": (xy[:, 1] / um_per_px).round(), "pxl_col_in_fullres": (xy[:, 0] / um_per_px).round()})
    (out / "spatial").mkdir(parents=True, exist_ok=True)
    positions.to_csv(out / "spatial" / "tissue_positions.csv", index=False)
    json.dump({"spot_diameter_fullres": 55 / um_per_px, "fiducial_diameter_fullres": 89 / um_per_px,
               "tissue_hires_scalef": 0.2, "tissue_lowres_scalef": 0.06}, open(out / "spatial" / "scalefactors_json.json", "w"))
    plt.imsave(out / "spatial" / "tissue_hires_image.png", np.ones((500, 500, 3)))
    plt.imsave(out / "spatial" / "tissue_lowres_image.png", np.ones((150, 150, 3)))

# %%
# 步骤：写 Visium HD 样本（16µm bin，每个 bin 1 个细胞，表达量较低）
# 输入：样本名、组别
# 输出：<样本>/outs/binned_outputs/square_016um/
def write_visium_hd(out, treated, size=800):
    bin_dir = out / "binned_outputs" / "square_016um"
    rows, cols = np.meshgrid(np.arange(size // 16), np.arange(size // 16), indexing="ij")
    xy = np.c_[cols.ravel() * 16.0, rows.ravel() * 16.0]
    counts = simulate_counts(layout(xy, treated, size)[:, None], 0.6)
    barcodes = [f"s_016um_{r:05d}_{c:05d}-1" for r, c in zip(rows.ravel(), cols.ravel())]
    write_10x_h5(bin_dir / "filtered_feature_bc_matrix.h5", counts, barcodes)
    um_per_px = 0.25
    (bin_dir / "spatial").mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"barcode": barcodes, "in_tissue": 1, "array_row": rows.ravel(), "array_col": cols.ravel(),
                  "pxl_row_in_fullres": xy[:, 1] / um_per_px, "pxl_col_in_fullres": xy[:, 0] / um_per_px}).to_parquet(bin_dir / "spatial" / "tissue_positions.parquet")
    json.dump({"spot_diameter_fullres": 16 / um_per_px, "bin_size_um": 16, "microns_per_pixel": um_per_px,
               "tissue_hires_scalef": 0.1, "tissue_lowres_scalef": 0.03}, open(bin_dir / "spatial" / "scalefactors_json.json", "w"))
    plt.imsave(bin_dir / "spatial" / "tissue_hires_image.png", np.ones((320, 320, 3)))

# %%
# 步骤：写 Xenium 样本（单细胞，坐标为 µm）
# 输入：样本名、组别
# 输出：<样本>/cell_feature_matrix.h5、cells.parquet
def write_xenium(out, treated, n=1500, size=1000):
    xy = rng.uniform(0, size, (n, 2))
    counts = simulate_counts(layout(xy, treated, size)[:, None], 1.0)
    cell_ids = [f"cell{i:05d}" for i in range(n)]
    write_10x_h5(out / "cell_feature_matrix.h5", counts, cell_ids)
    pd.DataFrame({"cell_id": cell_ids, "x_centroid": xy[:, 0], "y_centroid": xy[:, 1]}).to_parquet(out / "cells.parquet")

# %%
# 步骤：写 CosMx 样本（两个 FOV，全片像素坐标，含背景 cell_ID = 0 和阴性探针）
# 输入：样本名、组别
# 输出：<样本>/Run_exprMat_file.csv、Run_metadata_file.csv
def write_cosmx(out, treated, n=1500, size=1000):
    xy = rng.uniform(0, size, (n, 2))
    counts = simulate_counts(layout(xy, treated, size)[:, None], 1.0)
    fov = np.where(xy[:, 0] < size / 2, 1, 2)
    cell_id = np.concatenate([np.arange(1, (fov == 1).sum() + 1), np.arange(1, (fov == 2).sum() + 1)])
    order = np.argsort(fov, kind="stable")
    expr = pd.DataFrame(counts[order], columns=genes)
    expr.insert(0, "cell_ID", cell_id); expr.insert(0, "fov", fov[order])
    expr["NegPrb1"] = rng.poisson(0.1, n)
    background = pd.DataFrame([[1, 0] + [0] * (len(genes) + 1)], columns=expr.columns)
    out.mkdir(parents=True, exist_ok=True)
    pd.concat([expr, background]).to_csv(out / "Run_exprMat_file.csv", index=False)
    meta = expr[["fov", "cell_ID"]].copy()
    meta["CenterX_global_px"], meta["CenterY_global_px"] = xy[order, 0] / 0.12028, xy[order, 1] / 0.12028
    meta.to_csv(out / "Run_metadata_file.csv", index=False)

# %%
# 步骤：写 CODEX 样本（每行一个细胞的 marker 强度表，像素坐标）
# 输入：样本名、组别
# 输出：<样本>.csv
protein_markers = {"Tumor": ["PanCK", "Ki67"], "Fibroblasts": ["SMA", "Vimentin"], "Macrophage": ["CD68", "CD163"],
                   "T cells": ["CD3", "CD8"], "Endothelial": ["CD31", "CD34"]}
def write_codex(path, treated, n=2000, size=1000):
    xy = rng.uniform(0, size, (n, 2))
    labels = layout(xy, treated, size)
    channels = [m for ms in protein_markers.values() for m in ms] + ["DAPI"]
    intensity = rng.gamma(2, 20, (n, len(channels)))
    for i, cell_type in enumerate(cell_types):
        for m in protein_markers[cell_type]:
            intensity[labels == i, channels.index(m)] += 600
    table = pd.DataFrame(intensity.round(2), columns=channels)
    table.insert(0, "Centroid Y", xy[:, 1] / 0.377); table.insert(0, "Centroid X", xy[:, 0] / 0.377)
    path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(path, index=False)

# %%
# 步骤：写单细胞参考（每种细胞 400 个细胞，原始 counts）
# 输入：无
# 输出：<工作目录>/data/scrna_reference.h5ad
labels = np.repeat(np.arange(len(cell_types)), 400)
reference = sc.AnnData(sparse.csr_matrix(simulate_counts(labels[:, None], 1.5).astype(np.float32)), var=pd.DataFrame(index=genes))
reference.obs_names = [f"ref{i:05d}" for i in range(len(labels))]
reference.obs["cell_type"] = [cell_types[i] for i in labels]
(work_dir / "data").mkdir(parents=True, exist_ok=True)
reference.write_h5ad(work_dir / "data" / "scrna_reference.h5ad")

# %%
# 步骤：生成各平台的原始数据（低精度平台 2 个样本；单细胞级平台每组 2 个样本，供组间差异检验）
# 输入：上面的写入函数
# 输出：<工作目录>/data/<平台>/<样本>/
two_samples = [("S1", "Control", False), ("S2", "Treated", True)]
four_samples = [("S1", "Control", False), ("S2", "Control", False), ("S3", "Treated", True), ("S4", "Treated", True)]
sample_paths = {}
for platform, writer, sample_list in [("visium", write_visium, two_samples), ("visium_hd", write_visium_hd, four_samples),
                                      ("xenium", write_xenium, four_samples), ("cosmx", write_cosmx, four_samples)]:
    rows = []
    for sample, condition, treated in sample_list:
        sample_dir = work_dir / "data" / platform / sample / ("outs" if platform.startswith("visium") else "")
        writer(sample_dir, treated)
        rows.append({"sample": sample, "condition": condition, "path": str(sample_dir)})
    sample_paths[platform] = pd.DataFrame(rows)
rows = []
for sample, condition, treated in four_samples:
    path = work_dir / "data" / "codex" / f"{sample}.csv"
    write_codex(path, treated)
    rows.append({"sample": sample, "condition": condition, "path": str(path)})
sample_paths["codex"] = pd.DataFrame(rows)

# %%
# 步骤：建 7 个测试项目：配置来自 assets/configs/ 的预设，只改参考路径、线程数和 niche 数
# 输入：assets/configs/P*.yaml
# 输出：<工作目录>/projects/<流程>/ 下的 config.yaml、samples.tsv、markers.yaml
project_platform = {"P1_lowres_scrna": "visium", "P2_lowres_noref": "visium", "P3_hd_scrna": "visium_hd", "P4_hd_noref": "visium_hd",
                    "P5_imaging_scrna": "xenium", "P6_imaging_noref": "cosmx", "P7_protein": "codex"}
for pipeline, platform in project_platform.items():
    project = work_dir / "projects" / pipeline
    project.mkdir(parents=True, exist_ok=True)
    cfg = yaml.safe_load(open(skill_dir / "assets" / "configs" / f"{pipeline}.yaml"))
    cfg.update(platform=platform, threads=4, n_niches=4, scrna_reference=str(work_dir / "data" / "scrna_reference.h5ad") if cfg["scrna_reference"] else None)
    yaml.safe_dump(cfg, open(project / "config.yaml", "w"), allow_unicode=True, sort_keys=False)
    sample_paths[platform].to_csv(project / "samples.tsv", sep="\t", index=False)
    marker_table = protein_markers if platform == "codex" else {k: v[:3] for k, v in markers.items()}
    yaml.safe_dump(marker_table, open(project / "markers.yaml", "w"), allow_unicode=True)
print(work_dir)
