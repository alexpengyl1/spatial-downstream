# spatial-downstream

Downstream pipelines for spatial transcriptomics and spatial proteomics. For every step they name one method, say whether it runs in Python or R, and say whether samples can be analysed together or must be handled one at a time.

[中文版 README](README.zh-CN.md) · Version 1.0.0 · License: MIT · Languages: Python + R

Keywords: spatial transcriptomics, spatial proteomics, Visium, Visium HD, Stereo-seq, Slide-seq, Xenium, CosMx, MERSCOPE, CODEX, IMC, deconvolution, RCTD, cell2location, label transfer, spatial domain, BANKSY, CellCharter, cell niche, cellular neighbourhood, neighbourhood enrichment, MISTy, ligand–receptor, LIANA+, COMMOT, spatial CNV, scanpy, squidpy, FAIR.

## 1. Why this skill exists

Two things make downstream spatial analysis confusing.

The first is the number of tools. One task often has five candidates (RCTD or cell2location? BANKSY or CellCharter? squidpy, LIANA+, COMMOT or CellChat?), some exist only in R and some only in Python, and it is rarely clear which one suits a given platform or how to move data between the two languages.

The second is knowing which analyses can pool samples. Some steps can put every sample into one model; others only make sense inside one tissue section. Getting this wrong raises no error. The neighbour graph quietly links spots from different sections and the results are wrong.

This skill gives seven pipelines, P1 to P7, one for each combination of data type and whether a scRNA-seq reference is available; Figure 1 picks one from three questions. Each step uses one recommended method and states its language. Where a step needs R, a tested file-based bridge takes the data there and back. Figures 2 to 4 show what every step reads from the previous one and what it writes. Each step is also labelled as run jointly on all samples, run per sample, single-sample only, or needing several samples, with the reason (see "Multi-sample or single-sample?" below). The scripts are short, the environments are pinned, every run writes a provenance record, and all seven pipelines have been run end to end on simulated data in each platform's real file format.

Each step is a short top-to-bottom script split into `# %%` cells, so it runs in Jupyter or VS Code as well as with `python script.py`. Steps pass data through a documented AnnData contract, and `SKILL.md` tells the Claude Code or Codex assistant how to choose and adapt a pipeline.

## 2. Choose a pipeline

The documentation has four figures, from coarse to fine:

| Figure | File | Question it answers |
|---|---|---|
| Figure 1 · choose a pipeline | `docs/decision_en.svg` | Which of P1–P7 fits my data? |
| Figure 2 · step-by-step overview | `docs/steps_en.svg` | For each data type, what runs at every step? (all seven pipelines in one chart) |
| Figure 3 · single-pipeline diagram | `docs/P1_en.svg` … `docs/P7_en.svg` | How does the chosen pipeline run, script by script? What does each step read and write? |
| Figure 4 · data flow | `docs/dataflow_en.svg` | Which part of the previous file does each step use, and which fields does it add? (shared by all pipelines) |

Pick a pipeline in Figure 1, find its column in Figure 2, then open its Figure 3 and follow it. Figure 4 describes the files that connect the steps.

**Figure 1 · choose a pipeline** (answer the three questions from top to bottom):

![Figure 1: choose a pipeline](docs/decision_en.svg)

| ID | Data type | Platforms | scRNA-seq reference | Annotation | Languages |
|---|---|---|---|---|---|
| **P1** `P1_lowres_scrna` | spot-level (several cells per spot) | Visium, Stereo-seq bin ≥ 50 | yes | RCTD *full* (or cell2location) | Python + R (RCTD step) |
| **P2** `P2_lowres_noref` | spot-level | Visium, Stereo-seq bin ≥ 50 | no | marker scores | Python |
| **P3** `P3_hd_scrna` | near single-cell, sequencing-based | Visium HD 8/16 µm, Stereo-seq bin ≤ 20 / cellbin, Slide-seq | yes | RCTD *doublet* | Python + R (RCTD step) |
| **P4** `P4_hd_noref` | near single-cell, sequencing-based | as P3 | no | marker scores per cluster | Python |
| **P5** `P5_imaging_scrna` | single-cell, targeted gene panel | Xenium, CosMx, MERSCOPE | yes | Harmony + KNN label transfer | Python |
| **P6** `P6_imaging_noref` | single-cell, targeted gene panel | as P5 | no | marker scores per cluster | Python |
| **P7** `P7_protein` | single-cell, protein panel | CODEX/PhenoCycler, IMC, mIHC, MIBI | not applicable | marker scores + manual check | Python |

The three questions, in order:
1. Is one observation one cell? No → P1/P2. Yes → P3–P7.
2. Is it sequencing-based whole transcriptome (P3/P4), a targeted RNA panel (P5/P6) or protein (P7)?
3. Do you have a matched or same-tissue scRNA-seq reference with raw counts? Yes → odd-numbered pipeline; no → even-numbered.

Details and caveats for each case: [`references/decision_guide.md`](references/decision_guide.md).

## 3. What each pipeline runs

**Figure 2 · step-by-step overview of all seven pipelines.** Columns are data types and rows are steps in run order. Solid boxes run by default, dashed boxes are optional, grey boxes do not apply (the reason is written in the box), and a red "R" badge marks an R step.

![Figure 2: step-by-step overview](docs/steps_en.svg)

The same information as a table:

| Module | Step (script) | P1 | P2 | P3 | P4 | P5 | P6 | P7 | Language |
|---|---|---|---|---|---|---|---|---|---|
| Z | `validate_config.py` | ● | ● | ● | ● | ● | ● | ● | Py |
| A | `load_visium` / `load_visium_hd` / `load_stereo` / `load_xenium` / `load_cosmx` / `load_codex` | ● | ● | ● | ● | ● | ● | ● | Py |
| B | `preprocess_rna.py` (`preprocess_protein.py` for P7) | ● | ● | ● | ● | ● | ● | ● | Py |
| C | `prepare_reference` → `export_for_rctd` → **`run_rctd.Rmd`** → `import_rctd` | ● | | ● | | | | | Py → **R** → Py |
| C | `label_transfer.py` | | | | | ● | | | Py |
| C | `marker_annotation.py` | | ● | | ● | | ● | ● | Py |
| D | `domain_cellcharter.py` (or **`domain_banksy.Rmd`**) | ● | ● | ● | ● | ● | ● | ● | Py (or R) |
| D | `niche_composition.py` (niche / RCN / CN) | ● | ● | ● | ● | ● | ● | ● | Py |
| E | `misty_colocalization.py` | ● | ● | | | | | ● | Py |
| E | `nhood_enrichment.py` | | | ● | ● | ● | ● | ● | Py |
| E | `diff_nhood.py` (≥ 2 groups × ≥ 2 samples) | | | ● | ● | ● | ● | ● | Py |
| F | `lr_spatial_liana.py` | ● | ● | ● | ● | ● | ● | | Py |
| F | `ligrec_neighborhood.py`, `commot_flow.py` | opt | opt | opt | opt | opt | opt | | Py |
| G | `svg_moran.py` | ● | ● | ● | ● | ● | ● | ● | Py |
| G | `nmf_programs.py` | ● | ● | ● | ● | | | | Py |
| G | `cnv_infercnvpy.py`, `hotspot_getis.py`, `nmf_subtype.py`, `spatial_tendency_soapy.py` | opt | opt | opt | opt | opt | opt | opt | Py |
| Z | `record_provenance.py` | ● | ● | ● | ● | ● | ● | ● | Py |

● runs by default · opt = commented out in the pipeline; enable it after setting its parameters.

### Figure 3 · single-pipeline diagrams

Each diagram below is one column of Figure 2 drawn in detail. Under every step, "reads" is the part of the previous file the step uses and "writes" is the file it produces. Solid step badges run by default and dashed ones are optional. Blue marks Python steps, red marks R steps, yellow marks key outputs and the scRNA-seq reference, and grey is used for config, provenance and steps that do not apply. To regenerate the figures, run `python docs/make_figures.py`.

<details><summary><b>P1</b> · Low-resolution + scRNA-seq reference</summary>

![P1](docs/P1_en.svg)

</details>

<details><summary><b>P2</b> · Low-resolution, no reference</summary>

![P2](docs/P2_en.svg)

</details>

<details><summary><b>P3</b> · High-resolution sequencing + reference</summary>

![P3](docs/P3_en.svg)

</details>

<details><summary><b>P4</b> · High-resolution sequencing, no reference</summary>

![P4](docs/P4_en.svg)

</details>

<details><summary><b>P5</b> · Imaging single-cell + reference</summary>

![P5](docs/P5_en.svg)

</details>

<details><summary><b>P6</b> · Imaging single-cell, no reference</summary>

![P6](docs/P6_en.svg)

</details>

<details><summary><b>P7</b> · Spatial proteomics</summary>

![P7](docs/P7_en.svg)

</details>

### Figure 4 · data flow

Every step reads one part of the previous h5ad file and writes the next one. All downstream scripts read `05_niche.h5ad`. Field definitions are in [`schemas/anndata_contract.md`](schemas/anndata_contract.md).

![Figure 4: data flow](docs/dataflow_en.svg)

### Multi-sample or single-sample?

The rules below were checked against how the original example scripts handle samples. Steps that use only expression can pool samples, with Harmony correcting batches. Steps that use coordinates must keep neighbours within one tissue section. Statistics compared against a random baseline, such as permutation tests and Moran's I, are computed per sample and then summarised, because pooling mixes between-sample differences into them. In practice you put all samples in one `samples.tsv` and run the pipeline once. The full table, with what each original script does and why, is in [`references/multi_sample.md`](references/multi_sample.md).

| Handling | Steps | Why |
|---|---|---|
| Joint across samples (one model, labels comparable between samples) | QC & clustering (Harmony), cell2location, label transfer, marker annotation, CNV, spatial domain (CellCharter / BANKSY), niche / RCN / CN, LR within neighbourhoods, NMF programmes, micro-environment subtypes | Expression-only steps are pooled directly. Coordinate-based steps build neighbour graphs per sample and then cluster jointly, so niche and domain labels mean the same thing in every sample. The original RCN example does the same: neighbourhoods are counted per FOV and clustered together. |
| Per sample, then summarised (looped automatically; per-sample tables kept) | RCTD, neighbourhood enrichment (mean z and the number of samples that agree), MISTy (mean importance), spatially variable genes (mean Moran's I and the number of samples where it is significant), co-occurrence, LIANA+ spatial LR, hotspots and distances | The statistic is only meaningful within one section, or pooling would mix between-sample differences into it. In a test, a gene with only a batch shift got Moran's I = 0.19 when pooled and about 0 per sample. |
| Single sample (set `sample` at the top of the script and rerun) | COMMOT signal flow, spatial trends (SOAPy) | Optimal transport and region masks only make sense within one section. |
| Requires several samples | group differences in neighbourhood enrichment (`diff_nhood.py`); sample-level tests on `niche_composition_by_sample.csv` | Samples are the replicates, so you need at least 2 groups with 2 samples each; 3 per group is more reliable. |

## 4. Quick start

### Install the skill (Claude Code or Codex)

The skill is a folder with `SKILL.md` at the top, which is the format both Claude Code and OpenAI Codex read. Clone it from GitHub into the skills directory of the tool you use (this needs `git`):

```bash
# Claude Code (all projects)
git clone https://github.com/alexpengyl1/spatial-downstream.git ~/.claude/skills/spatial-downstream

# Codex
git clone https://github.com/alexpengyl1/spatial-downstream.git ~/.codex/skills/spatial-downstream
```

To use it in one Claude Code project only, clone it into `<project>/.claude/skills/spatial-downstream`. Without `git`, use Code → Download ZIP on the GitHub page, unzip it, rename the folder to `spatial-downstream` and move it into one of the skills directories above. To update later, run `git -C ~/.claude/skills/spatial-downstream pull`.

If you use both tools, install the skill once and link it into the other, e.g. `ln -s ~/.claude/skills/spatial-downstream ~/.codex/skills/spatial-downstream`. Restart the tool afterwards. It finds the skill through the `name` and `description` in `SKILL.md`, so you can simply ask, for example, "analyse my Visium data with a scRNA-seq reference".

The pipelines locate their scripts from their own path and run from either folder. You only need to set `SPATIAL_SKILL_DIR` when you run `validate_config` or `record_provenance` (script or notebook) on its own, outside a pipeline, and the skill is installed somewhere other than `~/.claude/skills/spatial-downstream`.

### Run a pipeline

```bash
SKILL_DIR=~/.claude/skills/spatial-downstream      # Codex: ~/.codex/skills/spatial-downstream

# 1. environments (once)
conda env create -f "$SKILL_DIR"/envs/spatial-py.yml
Rscript "$SKILL_DIR"/envs/install_r_packages.R      # only for P1 / P3 / BANKSY

# 2. project folder
mkdir my_project && cd my_project
cp "$SKILL_DIR"/assets/configs/P1_lowres_scrna.yaml config.yaml
cp "$SKILL_DIR"/assets/samples_template.tsv samples.tsv
cp "$SKILL_DIR"/assets/markers_template.yaml markers.yaml
#    edit config.yaml (species, reference path, radius…) and samples.tsv (sample, condition, path)

# 3. run
bash "$SKILL_DIR"/pipelines/P1_lowres_scrna.sh
```

`samples.tsv` has three tab-separated columns: `sample`, `condition` and `path`. The path is the Space Ranger `outs/` folder, the Xenium output folder, the CosMx folder or the CODEX CSV.

Each stage writes an h5ad file to `results/`, so you can restart from any stage by running the remaining scripts by hand.

### Example: four Visium samples with a scRNA-seq reference

Suppose you have four human tumour Visium samples (two control, two treated) and a scRNA-seq dataset from the same tissue in `ref/scrna.h5ad`, with raw counts and cell types in the `cell_type` column. You want to know whether the immune cells around the tumour change after treatment.

Visium is spot-level and you have a reference, so this is P1. You can ask the assistant exactly that question and give it the data paths, or set it up yourself:

```bash
SKILL_DIR=~/.claude/skills/spatial-downstream
mkdir tumor_visium && cd tumor_visium
cp "$SKILL_DIR"/assets/configs/P1_lowres_scrna.yaml config.yaml
cp "$SKILL_DIR"/assets/samples_template.tsv samples.tsv
cp "$SKILL_DIR"/assets/markers_template.yaml markers.yaml
```

Write one line per sample in `samples.tsv`:

```
sample	condition	path
C1	Control	/data/visium/C1/outs
C2	Control	/data/visium/C2/outs
T1	Treated	/data/visium/T1/outs
T2	Treated	/data/visium/T2/outs
```

In `config.yaml`, check `species: human`, `scrna_reference: ref/scrna.h5ad` and `scrna_celltype_key: cell_type`, and leave `n_niches` at 8 until you have seen the niche plots. Then run `bash "$SKILL_DIR"/pipelines/P1_lowres_scrna.sh` and look at `results/qc_violin.png` first; if the filtering thresholds do not suit your data, change them at the top of `preprocess_rna.py` and rerun from that step.

These files in `results/` answer the question:

| Question | Where to look |
|---|---|
| Which cells are in each spot? | `obsm['cell_props']` in `03_annotated.h5ad`; the RCTD report `reports/run_rctd.html` |
| Which micro-environments (niches) exist, and what are they made of? | `niche_composition_heatmap.png`, `spatial_niche.png` |
| How much of each sample does each niche cover, control vs treated? | `niche_by_sample.csv`, `niche_by_sample.png` |
| Does the same niche have a different make-up in different samples? | `niche_composition_by_sample.csv` |
| Which cell types sit next to each other? | `misty_juxta_mean.png` (mean over samples) |
| Ligand–receptor pairs between tumour and immune cells | `lr_local_top.png`, `lr_global_<sample>.csv` |

Niches are clustered on all four samples together, so niche 3 in C1 is the same niche as niche 3 in T2. RCTD, MISTy and Moran's I are computed per sample and the tables report the mean. With two samples per group, `niche_by_sample.csv` shows a trend but is too few for a reliable significance test; three or more samples per group is more reliable. Because spots are not single cells, P1 does not run neighbourhood enrichment; MISTy covers co-localisation instead.

### Interactive use: JupyterLab notebooks and R Markdown reports

Every Python script in `scripts/` has a JupyterLab notebook with the same code in `notebooks/`, in the same folder structure. A notebook opens with a description and a cell where you set `project_dir`. Each code cell starts with the three-line step / input / output header and is followed by a sentence saying what the step produced. The pipelines run the `.py` files, so after editing one, run `python tools/make_notebooks.py` to regenerate the notebooks and keep them identical.

The R steps are R Markdown reports (`run_rctd.Rmd`, `domain_banksy.Rmd`, `seurat_to_mtx.Rmd`) with numbered sections, a floating table of contents, code shown by default, an explanation before each chunk, check tables and `sessionInfo()` at the end. Knit them in RStudio after setting `project_dir` in the YAML `params`, or let the pipelines render them into `results/reports/`. Rendering needs pandoc, which ships with RStudio and Quarto; on the command line, set `RSTUDIO_PANDOC` if pandoc is not on the `PATH`.

## 5. Repository layout

```
spatial-downstream/
├── README.md, README.zh-CN.md   documentation (English / Chinese)
├── SKILL.md                     instructions for the Claude Code / Codex assistant
├── VERSION, CHANGELOG.md        version and history
├── CITATION.cff, codemeta.json  citation and machine-readable metadata
├── LICENSE                      MIT
├── pipelines/                   P1–P7 bash pipelines
├── notebooks/                   JupyterLab notebooks, one per Python script (generated)
├── tools/                       make_notebooks.py: regenerate notebooks from scripts
├── scripts/
│   ├── A_load/                  platform readers → 01_raw.h5ad
│   ├── B_preprocess/            QC, normalisation, Harmony, Leiden → 02_clustered.h5ad
│   ├── C_annotate/              RCTD, cell2location, label transfer, markers → 03_annotated.h5ad
│   ├── D_domain_niche/          spatial domains, niches → 04_domain.h5ad, 05_niche.h5ad
│   ├── E_neighborhood/          neighbourhood enrichment, group differences, MISTy
│   ├── F_communication/         ligand–receptor analyses
│   ├── G_special/               CNV, SVG, hotspots, NMF, subtypes, spatial trends
│   ├── bridge/                  R ↔ Python exchange
│   └── Z_provenance/            config validation, provenance record
├── docs/                        pipeline diagrams (SVG) and make_figures.py
├── assets/                      config presets (configs/P*.yaml), templates
├── schemas/                     config JSON Schema, AnnData data contract
├── envs/                        pinned Python env, R installer, tested versions
├── references/                  decision guide, method matrix, multi-sample guide, units, R–Python bridge, script map
└── tests/                       simulated-data generator and end-to-end test runner
```

## 6. Combining R and Python

The main object is always an h5ad file. When a step needs R (RCTD, BANKSY, Seurat-based tools), Python exports plain files, R returns a CSV keyed by barcode, and Python merges it back:

```
h5ad ──export (counts.mtx genes×cells, genes.txt, barcodes.txt, obs.csv, coords.csv)──▶ R ──result.csv──▶ merge by barcode ──▶ h5ad
```

Exchange raw integer counts only, use the same gene symbols on both sides, give coordinates in µm, and always merge by barcode rather than by row order. The full guide is [`references/r_python_bridge.md`](references/r_python_bridge.md).

## 7. Conventions

- Units: every radius and distance uses `obsm['spatial_um']` (µm). The readers convert platform units (pixels, DNB) to µm. See [`references/units_and_params.md`](references/units_and_params.md).
- Multiple samples: neighbour graphs never cross samples. They are built per sample (`library_key="sample"`), or the coordinates are offset per sample (MISTy, BANKSY), or the step loops over samples. See [`references/multi_sample.md`](references/multi_sample.md).
- Raw data: `layers['counts']` (RNA) and `layers['raw']` (protein) are never overwritten, and h5ad files keep all genes.
- Spot data: use cell-type proportions (niche, MISTy) rather than neighbourhood enrichment of the dominant label.
- Targeted panels: no CNV, pathway enrichment or RNA velocity. Check how many ligand–receptor pairs the panel covers before running LR analysis.
- macOS and Windows: Python's `spawn` start method makes process pools re-execute the whole script, so the scripts avoid them (progress bars are off in squidpy, and the CNV script sets `fork`).

## 8. Testing

```bash
bash "$SKILL_DIR"/tests/run_tests.sh /path/to/work_dir
```

`tests/make_test_data.py` writes small simulated datasets in each platform's real file format (Space Ranger HDF5 + spatial folder, Visium HD binned parquet, Xenium HDF5 + parquet, CosMx CSVs, CODEX CSV) plus a scRNA-seq reference. The runner then runs P1–P7 end to end and prints PASS or FAIL for each.

Tested: all seven pipelines on simulated data (macOS arm64, versions in [`envs/tested_versions.tsv`](envs/tested_versions.tsv)); readers for the Visium, Visium HD, Xenium, CosMx and CODEX formats; RCTD in full and doublet mode; the BANKSY round trip.
Not yet tested: real instrument data, `load_stereo.py` (stereopy), `run_cell2location.py`, `commot_flow.py`, `spatial_tendency_soapy.py`, Linux and Windows.

## 9. FAIR statement

| Principle | How this package addresses it |
|---|---|
| **F**indable | Descriptive name, keywords and version in `README`, `CITATION.cff` and `codemeta.json`; semantic versioning (`VERSION`, `CHANGELOG.md`). Archiving a tagged release on Zenodo gives a persistent identifier (DOI, see §10). |
| **A**ccessible | Plain-text scripts and documentation; MIT licence (free for any use, including commercial); only open-source dependencies from PyPI, CRAN and Bioconductor; no login or proprietary software. |
| **I**nteroperable | Community-standard formats only: AnnData `.h5ad`, Matrix Market `.mtx`, CSV/TSV, YAML, JSON Schema. A documented data contract ([`schemas/anndata_contract.md`](schemas/anndata_contract.md)) and a validated configuration ([`schemas/config.schema.json`](schemas/config.schema.json)); gene symbols (HGNC/MGI) and optional Cell Ontology IDs; coordinates in SI units (µm). |
| **R**eusable | Clear licence and citation; pinned, tested environments; every run writes `results/provenance/` (config copy, software versions, R session info, SHA-256 checksums of inputs and outputs); end-to-end tests on simulated data; each script states its inputs, outputs and parameters at the top. |

## 10. Citation and licence

Cite this package with the metadata in [`CITATION.cff`](CITATION.cff), and cite the methods you used:

| Step | Method |
|---|---|
| preprocessing | Scanpy (Wolf et al., *Genome Biol* 2018); Harmony (Korsunsky et al., *Nat Methods* 2019) |
| spatial statistics | Squidpy (Palla et al., *Nat Methods* 2022) |
| deconvolution | RCTD (Cable et al., *Nat Biotechnol* 2022); cell2location (Kleshchevnikov et al., *Nat Biotechnol* 2022) |
| spatial domains | BANKSY (Singhal et al., *Nat Genet* 2024); CellCharter (Varrone et al., *Nat Genet* 2024) |
| co-localisation, ligand–receptor | MISTy (Tanevski et al., *Genome Biol* 2022); LIANA+ (Dimitrov et al., *Nat Cell Biol* 2024); CellPhoneDB (Efremova et al., *Nat Protoc* 2020); COMMOT (Cang et al., *Nat Methods* 2023) |
| CNV | inferCNV (Tirosh et al., *Science* 2016) / infercnvpy |
| hotspots | Getis & Ord (1992), PySAL `esda` |
| spatial trends | SOAPy |

To get a DOI, enable the Zenodo integration for the GitHub repository and publish a release tagged `v1.0.0`, then add the DOI to `CITATION.cff` and `codemeta.json`.

Licence: [MIT](LICENSE). Copyright © 2026 Yueling Peng. Anyone may use, copy, modify and redistribute this package for any purpose, including commercial use, as long as the copyright and licence notice is kept.
