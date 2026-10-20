# AnnData data contract / AnnData 数据约定

Every script reads and writes AnnData (`.h5ad`) objects that follow this contract, so steps can be swapped or reused in other pipelines.
所有脚本读写的 h5ad 都遵守下面的约定，因此步骤之间可以替换，也可以在别的流程里复用。

## Stage files / 阶段文件

Diagram / 图示：[`docs/dataflow_en.svg`](../docs/dataflow_en.svg) · [`docs/dataflow_zh.svg`](../docs/dataflow_zh.svg)

| file | written by | adds |
|---|---|---|
| `01_raw.h5ad` | `A_load/*` | raw data, coordinates |
| `02_clustered.h5ad` | `B_preprocess/*` | normalised X, PCA, Leiden |
| `03_annotated.h5ad` | `C_annotate/*` | cell type / composition |
| `04_domain.h5ad` | `D_domain_niche/domain_*` (+ `bridge/import_obs_csv.py`) | spatial domain |
| `05_niche.h5ad` | `D_domain_niche/niche_composition.py` | niche, neighbourhood composition |
| tables & figures | `E_*`, `F_*`, `G_*` | results only, no new h5ad |

## Keys / 字段

| slot | key | type | meaning | required from |
|---|---|---|---|---|
| `obs` | `sample` | category | sample id (matches `samples.tsv`) | 01 |
| `obs` | `condition` | category | experimental group | 01 |
| `obsm` | `spatial` | float (n × 2) | platform-native coordinates (pixels, DNB, µm) | 01 |
| `obsm` | `spatial_um` | float (n × 2) | coordinates in **micrometres**; all radii use this | 01 |
| `layers` | `counts` | sparse int-valued | raw counts (RNA platforms) | 01 |
| `layers` | `raw` | float | raw intensities (protein platforms) | 01 |
| `X` | — | float | log1p-normalised (RNA) / arcsinh z-score (protein) | 02 |
| `var` | `highly_variable` | bool | genes used for PCA (all genes for targeted panels) | 02 |
| `obsm` | `X_pca`, `X_pca_harmony` | float | embeddings (harmony when > 1 sample) | 02 |
| `obs` | `leiden` | category | expression clusters | 02 |
| `obs` | `cell_type` | category | cell type (dominant type for spot data) | 03 |
| `obs` | `cell_type_score` | float | label-transfer confidence | 03 (label transfer) |
| `obsm` | `cell_props` | DataFrame (n × types), rows sum to 1 | deconvolution proportions | 03 (RCTD / cell2location) |
| `obsm` | `cell_scores` | DataFrame (n × types) | marker scores (not proportions) | 03 (markers) |
| `obs` | `domain` | category | spatial domain | 04 |
| `obs` | `niche` | category | niche / RCN / CN | 05 |
| `obsm` | `nhood_composition` | DataFrame (n × types), rows sum to 1 | neighbourhood cell-type composition | 05 |

## Exchange files for R / 与 R 交换的文件

`bridge/h5ad_to_mtx.py` and `C_annotate/export_for_rctd.py` write:

| file | format |
|---|---|
| `counts.mtx` | Matrix Market, **genes × cells**, raw integer-valued counts |
| `genes.txt` | one gene symbol per line |
| `barcodes.txt` | one cell / spot id per line (unique across samples) |
| `obs.csv` | metadata, first column = barcode |
| `coords.csv` | `x`, `y` in µm, first column = barcode |

R scripts return results as CSV with the barcode in the first column; Python merges them back by barcode (never by row order).
R 脚本以 CSV 返回结果，第一列为条形码；Python 按条形码合并（不按行号）。

## Naming / 命名

- Gene identifiers: HGNC (human) or MGI (mouse) symbols in `var_names`; keep Ensembl IDs in `var['gene_ids']` when available.
- Cell type names: free text without `/` (RCTD restriction). For cross-study reuse, add an optional `obs['cell_type_ontology_id']` with Cell Ontology IDs (e.g. `CL:0000235` macrophage).
- 基因名用 symbol；细胞类型名里不要有 `/`；跨研究复用时建议加一列 Cell Ontology ID。
