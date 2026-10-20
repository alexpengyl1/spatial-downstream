# 生成流程图（SVG，中英文各一份）/ Generate pipeline figures (SVG, English and Chinese)
# 运行 / usage：python docs/make_figures.py   → 图 1 decision_*、图 2 steps_*、图 3 P1…P7_*、图 4 dataflow_*（en / zh，SVG）
# 风格参考 pysc 的流程图；配色为低饱和度的 BioRender 风格。只用标准库，不需要额外安装。

# %%
# 步骤：导入、输出目录、配色和字体
# 输入：无
# 输出：OUT、PALETTE、FONT
from pathlib import Path
from html import escape

OUT = Path(__file__).resolve().parent
FONT = "-apple-system,'PingFang SC','Noto Sans CJK SC','Segoe UI','Helvetica Neue',Arial,sans-serif"
# 四种颜色：蓝、红、黄 + 中性灰（低饱和度 BioRender 风格）；每种颜色只取同一色相的深浅几档
BLUE = ("#5E81AC", "#EAF0F7")                  # 蓝：Python 分析步骤
RED = ("#C47A72", "#F8EAE8")                   # 红：R 步骤（RCTD、BANKSY）
YELLOW = ("#C9A24A", "#FAF3DE")                # 黄：关键输出文件、"有单细胞参考"
YELLOW_DARK = "#86692A"                        # 黄深一档：黄底上的文字
GREY = ("#7D8A99", "#EEF1F4")                  # 灰：输入、配置、溯源、不适用、连线
GREY_DARK, GREY_LINE = "#3B4252", "#D8DEE6"    # 灰的深浅两档：正文文字、边框
TEXT, SUB, BORDER, CANVAS = GREY_DARK, GREY[0], GREY_LINE, "#FFFFFF"
PALETTE = {
    "input": GREY, "Z": GREY, "R": RED,
    "A": BLUE, "B": BLUE, "C": BLUE, "D": BLUE, "E": BLUE, "F": BLUE, "G": BLUE, "Py": BLUE,
    "acc": YELLOW,
}

# %%
# 步骤：基础绘图函数（文字、矩形、标签、箭头、图标）
# 输入：坐标和样式
# 输出：SVG 字符串片段
def width_of(s, size):
    return sum(size * (1.0 if ord(c) > 0x2E80 else 0.58) for c in s)

def text(x, y, s, size=14, weight=400, fill=TEXT, anchor="start"):
    return f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" font-weight="{weight}" fill="{fill}" text-anchor="{anchor}">{escape(s)}</text>'

def rect(x, y, w, h, fill=CANVAS, stroke=BORDER, rx=16, sw=1.5, dash=None, shadow=False):
    extra = f' stroke-dasharray="{dash}"' if dash else ""
    filt = ' filter="url(#shadow)"' if shadow else ""
    return f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" rx="{rx}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"{extra}{filt}/>'

def pill(x, y, label, key, size=12, solid=False, dashed=False):
    color, light = PALETTE[key]
    w = width_of(label, size) + 18
    fill, fg = (color, "#FFFFFF") if solid else (light, color)
    dash = ' stroke-dasharray="3 3"' if dashed else ""
    svg = (f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{size + 10}" rx="{(size + 10) / 2}" fill="{fill}" stroke="{color}" stroke-width="1"{dash}/>'
           + text(x + w / 2, y + size + 3, label, size, 600, fg, "middle"))
    return svg, w

def chevron(x, y):
    return f'<path d="M{x - 5},{y - 9} L{x + 4},{y} L{x - 5},{y + 9}" fill="none" stroke="{GREY[0]}" stroke-width="2.4" stroke-linecap="round" stroke-linejoin="round"/>'

def arrow_down(x, y1, y2, color=GREY[0]):
    return (f'<line x1="{x}" y1="{y1}" x2="{x}" y2="{y2 - 9}" stroke="{color}" stroke-width="2"/>'
            f'<path d="M{x - 7},{y2 - 11} L{x},{y2} L{x + 7},{y2 - 11} Z" fill="{color}"/>')

def icon(kind, cx, cy, key, r=26):
    color, light = PALETTE[key]
    s = f'<circle cx="{cx}" cy="{cy}" r="{r}" fill="{light}" stroke="{color}" stroke-width="1.8"/>'
    g = f'<g fill="none" stroke="{color}" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" transform="translate({cx},{cy})">'
    glyph = {
        "data": '<ellipse cx="0" cy="-7" rx="10" ry="4"/><path d="M-10,-7 v14 a10,4 0 0 0 20,0 v-14"/><path d="M-10,0 a10,4 0 0 0 20,0"/>',
        "check": '<rect x="-10" y="-11" width="20" height="22" rx="3"/><path d="M-5,0 l3.5,3.5 l6.5,-7"/>',
        "load": '<path d="M-8,-11 h10 l6,6 v16 h-16 z"/><path d="M2,-11 v6 h6"/><path d="M-4,3 h8 M-4,7 h6"/>',
        "cluster": '<circle cx="-6" cy="-5" r="2.5" fill="currentColor"/><circle cx="6" cy="-6" r="2.5"/><circle cx="0" cy="6" r="2.5"/><circle cx="-9" cy="5" r="2"/><circle cx="9" cy="5" r="2"/>',
        "tag": '<path d="M-10,-3 v-8 h8 l12,12 -8,8 z"/><circle cx="-5" cy="-6" r="1.8"/>',
        "bridge": '<path d="M-10,-4 h18 l-4,-4 M10,4 h-18 l4,4"/>',
        "domain": '<rect x="-11" y="-11" width="10" height="10" rx="2"/><rect x="1" y="-11" width="10" height="10" rx="2"/><rect x="-11" y="1" width="10" height="10" rx="2"/><rect x="1" y="1" width="10" height="10" rx="2"/>',
        "niche": '<path d="M0,-12 l10.4,6 v12 l-10.4,6 l-10.4,-6 v-12 z"/><circle cx="0" cy="0" r="3"/>',
        "nhood": '<circle cx="0" cy="0" r="3.5"/><circle cx="-10" cy="-8" r="2.5"/><circle cx="10" cy="-8" r="2.5"/><circle cx="-10" cy="8" r="2.5"/><circle cx="10" cy="8" r="2.5"/><path d="M-3,-2 L-8,-6 M3,-2 L8,-6 M-3,2 L-8,6 M3,2 L8,6"/>',
        "comm": '<circle cx="-8" cy="0" r="5"/><circle cx="9" cy="0" r="5"/><path d="M-2,0 h7 l-3,-3 M5,0 l-3,3"/>',
        "special": '<path d="M-10,9 v-6 M-4,9 v-12 M2,9 v-8 M8,9 v-16"/>',
        "prov": '<path d="M0,-12 l10,4 v7 c0,7 -5,11 -10,13 c-5,-2 -10,-6 -10,-13 v-7 z"/><path d="M-4,0 l3,3 l5,-6"/>',
        "file": '<path d="M-9,-12 h11 l7,7 v17 h-18 z"/><path d="M2,-12 v7 h7"/>',
    }[kind]
    return s + g.replace("currentColor", color) + glyph + "</g>"

# %%
# 步骤：界面文字（中英文）
# 输入：无
# 输出：UI
UI = {
    "en": dict(upstream="Upstream · load, QC, annotation", structure="Spatial structure", downstream="Downstream · all read 05_niche.h5ad",
               optional="optional", na="Not applicable", prov="Provenance", prov_sub="record_provenance.py → results/provenance/: config copy, software versions, SHA-256 of inputs and outputs",
               niche_out="input to all downstream", ref_yes="scRNA-seq reference: yes", ref_no="scRNA-seq reference: no", lang="Languages",
               E="Neighbourhoods", F="Communication", G="Special analyses", detail="Annotation detail", alt="alternative",
               spans="reference not applicable"),
    "zh": dict(upstream="上游 · 读入、质控、注释", structure="空间结构", downstream="下游 · 都读取 05_niche.h5ad",
               optional="可选", na="不适用", prov="溯源记录", prov_sub="record_provenance.py → results/provenance/：配置副本、软件版本、输入和输出文件的 SHA-256",
               niche_out="下游分析的共同输入", ref_yes="单细胞参考：有", ref_no="单细胞参考：无", lang="语言",
               E="邻域与共定位", F="细胞通讯", G="个性化分析", detail="注释细节", alt="替代方案",
               spans="单细胞参考不适用"),
}

# %%
# 步骤：7 套流程的内容（每个步骤：英文名、中文名、脚本、语言、是否可选）
# 输入：无
# 输出：PIPELINES
def step(en, zh, script, lang="Py", optional=False, out=""):
    return dict(en=en, zh=zh, script=script, lang=lang, optional=optional, out=out)

# 下游脚本：从 05_niche.h5ad 读哪一部分（英文、中文），写出什么
SCRIPT_IO = {
    "misty_colocalization.py": ("composition + spatial_um, per sample", "细胞组成 + spatial_um（逐样本）", "misty_juxta_mean.csv"),
    "nhood_enrichment.py": ("cell_type + spatial_um, per sample", "cell_type + spatial_um（逐样本）", "nhood_enrichment_mean.csv"),
    "diff_nhood.py": ("cell_type + condition + spatial_um", "cell_type + 分组 + spatial_um", "diff_nhood_*.csv / .png"),
    "lr_spatial_liana.py": ("X (log) + spatial_um", "X（log 表达）+ spatial_um", "lr_global_*.csv"),
    "commot_flow.py": ("X (log) + spatial_um, one sample", "X + spatial_um（单样本）", "commot_*.png"),
    "ligrec_neighborhood.py": ("X + cell_type + spatial_um", "X + cell_type + spatial_um", "ligrec_<type>.csv"),
    "svg_moran.py": ("X + spatial_um, per sample", "X + spatial_um（逐样本）", "svg_moran.csv"),
    "nmf_programs.py": ("X (HVGs) + niche", "X（高变基因）+ niche", "nmf_top_genes.csv"),
    "cnv_infercnvpy.py": ("X + cell_type + GTF", "X + cell_type + GTF", "cnv_per_spot.csv"),
    "hotspot_getis.py": ("composition (obsm) + spatial_um", "obsm 细胞组成 + spatial_um", "hotspot_distance.csv"),
    "spatial_tendency_soapy.py": ("X + niche, one sample", "X + niche（单样本）", "tendency_*.csv"),
    "nmf_subtype.py": ("obsm nhood_composition + cell_type", "obsm nhood_composition + cell_type", "subtype_markers.csv"),
}

RCTD_DETAIL = [step("Prepare reference", "准备参考", "prepare_reference.py", out="ref_prepared.h5ad"),
               step("Export files", "导出文件", "export_for_rctd.py", out="rctd/*.mtx, coords.csv"),
               step("RCTD", "RCTD", "run_rctd.Rmd", "R", out="rctd/*/weights.csv"),
               step("Import results", "读回结果", "import_rctd.py", out="03_annotated.h5ad")]
LOWRES_EFG = dict(
    E=[step("MISTy co-localisation", "MISTy 共定位", "misty_colocalization.py")],
    F=[step("Spatial ligand–receptor", "空间配受体", "lr_spatial_liana.py"), step("Signal flow (COMMOT)", "信号流向（COMMOT）", "commot_flow.py", optional=True)],
    G=[step("Spatially variable genes", "空间可变基因", "svg_moran.py"), step("NMF programmes", "NMF 表达程序", "nmf_programs.py"),
       step("Spatial CNV", "空间 CNV", "cnv_infercnvpy.py", optional=True), step("Hotspots & distances", "热点与距离", "hotspot_getis.py", optional=True),
       step("Spatial trends", "空间趋势", "spatial_tendency_soapy.py", optional=True)],
    na=("neighbourhood enrichment of dominant labels", "用主导类型做邻域富集"))
HD_EFG = dict(
    E=[step("Neighbourhood enrichment", "邻域富集", "nhood_enrichment.py"), step("Group differences", "组间差异", "diff_nhood.py")],
    F=[step("Spatial ligand–receptor", "空间配受体", "lr_spatial_liana.py"), step("LR within neighbourhoods", "邻域内配受体", "ligrec_neighborhood.py", optional=True)],
    G=[step("Spatially variable genes", "空间可变基因", "svg_moran.py"), step("NMF programmes", "NMF 表达程序", "nmf_programs.py"),
       step("Micro-environment subtypes", "微环境亚型", "nmf_subtype.py", optional=True), step("Spatial CNV", "空间 CNV", "cnv_infercnvpy.py", optional=True)],
    na=None)
IMAGING_EFG = dict(
    E=[step("Neighbourhood enrichment", "邻域富集", "nhood_enrichment.py"), step("Group differences", "组间差异", "diff_nhood.py")],
    F=[step("Spatial ligand–receptor", "空间配受体", "lr_spatial_liana.py"), step("LR within neighbourhoods", "邻域内配受体", "ligrec_neighborhood.py", optional=True)],
    G=[step("Spatially variable genes", "空间可变基因", "svg_moran.py"), step("Micro-environment subtypes", "微环境亚型", "nmf_subtype.py", optional=True),
       step("Hotspots & distances", "热点与距离", "hotspot_getis.py", optional=True)],
    na=("CNV, pathway enrichment, COMMOT", "CNV、通路富集、COMMOT"))
PROTEIN_EFG = dict(
    E=[step("Neighbourhood enrichment", "邻域富集", "nhood_enrichment.py"), step("Group differences", "组间差异", "diff_nhood.py"),
       step("MISTy co-localisation", "MISTy 共定位", "misty_colocalization.py")],
    F=[],
    G=[step("Spatially variable markers", "空间可变 marker", "svg_moran.py"), step("Micro-environment subtypes", "微环境亚型", "nmf_subtype.py", optional=True),
       step("Hotspots & distances", "热点与距离", "hotspot_getis.py", optional=True)],
    na=("deconvolution, CNV, ligand–receptor", "去卷积、CNV、配受体"))

def pipeline(pid, en, zh, platforms_en, platforms_zh, ref, raw_en, raw_zh, loader, pre, annot, detail, niche_en, niche_zh, efg, alt=None):
    return dict(id=pid, en=en, zh=zh, platforms_en=platforms_en, platforms_zh=platforms_zh, ref=ref, raw_en=raw_en, raw_zh=raw_zh,
                loader=loader, pre=pre, annot=annot, detail=detail, niche_en=niche_en, niche_zh=niche_zh, efg=efg, alt=alt)

PIPELINES = [
    pipeline("P1", "Low-resolution spatial + scRNA-seq reference", "低精度空转 + 单细胞参考", "Visium · Stereo-seq bin ≥ 50", "Visium · Stereo-seq bin ≥ 50", True,
             "Space Ranger outs/ · GEF", "Space Ranger outs/ · GEF", "load_visium.py", "preprocess_rna.py",
             step("Deconvolution", "去卷积", "RCTD full mode", "R"), RCTD_DETAIL, "Cell niche", "细胞 niche", LOWRES_EFG,
             alt=step("cell2location (GPU)", "cell2location（GPU）", "run_cell2location.py")),
    pipeline("P2", "Low-resolution spatial, no reference", "低精度空转，无单细胞参考", "Visium · Stereo-seq bin ≥ 50", "Visium · Stereo-seq bin ≥ 50", False,
             "Space Ranger outs/ · GEF", "Space Ranger outs/ · GEF", "load_visium.py", "preprocess_rna.py",
             step("Marker scores", "marker 打分", "marker_annotation.py"), None, "Niche from scores", "基于打分的 niche", LOWRES_EFG),
    pipeline("P3", "High-resolution sequencing + scRNA-seq reference", "高精度测序型 + 单细胞参考", "Visium HD · Stereo-seq bin ≤ 20 · Slide-seq", "Visium HD · Stereo-seq bin ≤ 20 · Slide-seq", True,
             "binned_outputs/square_016um", "binned_outputs/square_016um", "load_visium_hd.py", "preprocess_rna.py",
             step("RCTD annotation", "RCTD 注释", "RCTD doublet mode", "R"), RCTD_DETAIL, "Niche / RCN", "niche / RCN", HD_EFG,
             alt=step("Label transfer", "标签转移（分割细胞）", "label_transfer.py")),
    pipeline("P4", "High-resolution sequencing, no reference", "高精度测序型，无单细胞参考", "Visium HD · Stereo-seq bin ≤ 20 · Slide-seq", "Visium HD · Stereo-seq bin ≤ 20 · Slide-seq", False,
             "binned_outputs/square_016um", "binned_outputs/square_016um", "load_visium_hd.py", "preprocess_rna.py",
             step("Cluster annotation", "按聚类注释", "marker_annotation.py"), None, "Niche / RCN", "niche / RCN", HD_EFG),
    pipeline("P5", "Imaging single-cell + scRNA-seq reference", "成像型单细胞 + 单细胞参考", "Xenium · CosMx · MERSCOPE", "Xenium · CosMx · MERSCOPE", True,
             "Xenium / CosMx output", "Xenium / CosMx 输出", "load_xenium.py", "preprocess_rna.py",
             step("Label transfer", "标签转移", "Harmony + KNN"),
             [step("Prepare reference", "准备参考", "prepare_reference.py", out="ref_prepared.h5ad"), step("Harmony + KNN", "Harmony + KNN", "label_transfer.py", out="03_annotated.h5ad")],
             "RCN", "RCN", IMAGING_EFG),
    pipeline("P6", "Imaging single-cell, no reference", "成像型单细胞，无单细胞参考", "Xenium · CosMx · MERSCOPE", "Xenium · CosMx · MERSCOPE", False,
             "Xenium / CosMx output", "Xenium / CosMx 输出", "load_xenium.py", "preprocess_rna.py",
             step("Cluster annotation", "按聚类注释", "marker_annotation.py"), None, "RCN", "RCN", IMAGING_EFG),
    pipeline("P7", "Spatial proteomics", "空间蛋白组", "CODEX · IMC · mIHC · MIBI", "CODEX · IMC · mIHC · MIBI", None,
             "single-cell intensity CSV", "单细胞强度 CSV", "load_codex.py", "preprocess_protein.py",
             step("Marker annotation", "marker 注释", "marker_annotation.py"), None, "CN", "CN", PROTEIN_EFG),
]

# %%
# 步骤：画一个上游节点（图标 + 标题 + 脚本名 + 语言标签）
# 输入：中心坐标、步骤、模块、语言
# 输出：SVG 片段
def node(cx, y, kind, key, title, sub, lang_key=None, io=None):
    s = icon(kind, cx, y + 30, key)
    s += text(cx, y + 82, title, 15, 600, TEXT, "middle") + text(cx, y + 101, sub, 12, 400, SUB, "middle")
    if lang_key:
        label = "R" if lang_key == "R" else "Python"
        w = width_of(label, 11) + 18
        p, _ = pill(cx - w / 2, y + 110, label, lang_key, 11)
        s += p
    if io:                                   # io = (读入, 写出)：读上一步的哪部分数据、写出什么文件
        s += text(cx, y + 150, io[0], 10.5, 400, SUB, "middle") + text(cx, y + 166, io[1], 10.5, 600, YELLOW_DARK, "middle")
    return s

def upstream_io(p, lang):
    zh = lang == "zh"
    r, w = ("读 ", "写 ") if zh else ("reads ", "writes ")
    protein = p["id"] == "P7"
    annot_in = {"R": "02_clustered · counts + 参考" if zh else "02_clustered · counts + ref"}.get(p["annot"]["lang"],
               ("02_clustered · counts + 参考" if zh else "02_clustered · counts + ref") if p["id"] == "P5" else ("02_clustered · X + markers"))
    niche_in = {"P1": "04_domain · obsm cell_props", "P2": "04_domain · obsm cell_scores"}.get(p["id"], "04_domain · cell_type + 坐标" if zh else "04_domain · cell_type + coords")
    return dict(
        raw=("—", w + ("仪器原始输出" if zh else "instrument output")),
        check=(r + "config.yaml · samples.tsv", w + ("无（只检查）" if zh else "nothing (check only)")),
        load=(r + ("原始输出目录" if zh else "raw output folder"), w + "01_raw.h5ad"),
        qc=(r + ("01_raw · layers raw" if protein else "01_raw · layers counts"), w + "02_clustered.h5ad"),
        annot=(r + annot_in, w + "03_annotated.h5ad"),
        domain=(r + (("03_annotated · X + 坐标" if protein else "03_annotated · X_pca + 坐标") if zh else ("03_annotated · X + coords" if protein else "03_annotated · X_pca + coords")), w + "04_domain.h5ad"),
        niche=(r + niche_in, w + "05_niche.h5ad"))

# %%
# 步骤：画一张流程图
# 输入：一套流程的内容、语言
# 输出：SVG 字符串
def draw_pipeline(p, lang):
    u = UI[lang]
    W, M = 1100, 30
    parts, y = [], 0
    title = f"图 3 · {p['id']} · {p[lang]}" if lang == "zh" else f"Figure 3 · {p['id']} · {p[lang]}"
    parts.append(text(M + 6, 52, title, 26, 700))
    lane = {"P1": ("spot 级", "Spot-level"), "P2": ("spot 级", "Spot-level"), "P3": ("近单细胞，测序型", "High-res sequencing"), "P4": ("近单细胞，测序型", "High-res sequencing"),
            "P5": ("单细胞，靶向 RNA", "Imaging, targeted RNA"), "P6": ("单细胞，靶向 RNA", "Imaging, targeted RNA"), "P7": ("空间蛋白组", "Spatial proteomics")}[p["id"]]
    where = (f"{p['platforms_zh']} · 对应图 2 中「{lane[0]}」一列" if lang == "zh" else f"{p['platforms_en']} · the “{lane[1]}” column in Figure 2")
    parts.append(text(M + 6, 80, where, 14, 400, SUB))
    x = M + 6
    for label, key in [(u["ref_yes"] if p["ref"] else (u["ref_no"] if p["ref"] is False else u["spans"]), "C"),
                       (f"{u['lang']}: Python + R" if p["annot"]["lang"] == "R" else f"{u['lang']}: Python", "R" if p["annot"]["lang"] == "R" else "Py")]:
        svg, w = pill(x, 94, label, key, 12)
        parts.append(svg); x += w + 10
    y = 136

    # 上游卡片
    has_detail = p["detail"] is not None
    io = upstream_io(p, lang)
    card_h = 236 + (134 if has_detail else 0)
    parts.append(rect(M, y, W - 2 * M, card_h, shadow=True))
    svg, _ = pill(M + 22, y + 16, u["upstream"], "A", 13)
    parts.append(svg)
    nodes = [("data", "input", "Raw data" if lang == "en" else "原始数据", p["raw_" + lang], None, io["raw"]),
             ("check", "Z", "Check config" if lang == "en" else "检查配置", "validate_config.py", "Py", io["check"]),
             ("load", "A", "Load" if lang == "en" else "读入", p["loader"], "Py", io["load"]),
             ("cluster", "B", "QC & clustering" if lang == "en" else "质控与聚类", p["pre"], "Py", io["qc"]),
             ("tag", "C", p["annot"][lang], p["annot"]["script"], p["annot"]["lang"], io["annot"])]
    xs = [M + 110 + i * 210 for i in range(5)]
    for i, (kind, key, t, s, lk, node_io) in enumerate(nodes):
        parts.append(node(xs[i], y + 46, kind, key, t, s, lk, node_io))
        if i < 4:
            parts.append(chevron((xs[i] + xs[i + 1]) / 2, y + 76))
    if has_detail:
        dy = y + 242
        dx0, dx1 = M + 22, W - M - 22
        parts.append(rect(dx0, dy, dx1 - dx0, 104, fill=PALETTE["R"][1] if p["annot"]["lang"] == "R" else PALETTE["C"][1],
                          stroke=PALETTE["R"][0] if p["annot"]["lang"] == "R" else PALETTE["C"][0], rx=12, sw=1.2, dash="5 4"))
        parts.append(text(dx0 + 16, dy + 24, u["detail"], 13, 600, PALETTE["R"][0] if p["annot"]["lang"] == "R" else PALETTE["C"][0]))
        n = len(p["detail"])
        span = (dx1 - dx0 - 40) / (n + (1 if p["alt"] else 0))
        for i, st in enumerate(p["detail"]):
            cx = dx0 + 20 + span * (i + 0.5)
            parts.append(text(cx, dy + 46, st[lang], 14, 600, TEXT, "middle"))
            svg, _ = pill(cx + width_of(st[lang], 14) / 2 + 6, dy + 32, "R" if st["lang"] == "R" else "Py", "R" if st["lang"] == "R" else "Py", 10, solid=st["lang"] == "R")
            parts.append(svg + text(cx, dy + 66, st["script"], 12, 400, SUB, "middle"))
            parts.append(text(cx, dy + 86, ("写 " if lang == "zh" else "writes ") + st["out"], 11, 600, YELLOW_DARK, "middle"))
            if i < n - 1:
                parts.append(chevron(cx + span / 2, dy + 62))
        if p["alt"]:
            cx = dx0 + 20 + span * (n + 0.5)
            parts.append(f'<line x1="{cx - span / 2}" y1="{dy + 30}" x2="{cx - span / 2}" y2="{dy + 70}" stroke="{BORDER}" stroke-width="1.5"/>')
            parts.append(text(cx, dy + 26, u["alt"], 11, 600, SUB, "middle"))
            parts.append(text(cx, dy + 46, p["alt"][lang], 13, 600, TEXT, "middle"))
            parts.append(text(cx, dy + 66, p["alt"]["script"], 12, 400, SUB, "middle"))
    y += card_h
    parts.append(arrow_down(W / 2, y + 2, y + 30))
    y += 34

    # 空间结构卡片
    parts.append(rect(M, y, W - 2 * M, 252, shadow=True))
    svg, _ = pill(M + 22, y + 16, u["structure"], "D", 13)
    parts.append(svg)
    parts.append(node(M + 110, y + 46, "domain", "D", "Spatial domain" if lang == "en" else "空间域", "domain_cellcharter.py", "Py", io["domain"]))
    parts.append(text(M + 110, y + 236, ("or domain_banksy.Rmd" if lang == "en" else "或 domain_banksy.Rmd"), 12, 600, PALETTE["R"][0], "middle"))
    parts.append(chevron(M + 215, y + 76))
    parts.append(node(M + 320, y + 46, "niche", "D", p["niche_" + lang], "niche_composition.py", "Py", io["niche"]))
    parts.append(chevron(M + 425, y + 76))
    bx, bw = M + 470, W - 2 * M - 492
    parts.append(rect(bx, y + 64, bw, 110, fill=PALETTE["acc"][1], stroke=PALETTE["acc"][0], rx=16, sw=2))
    parts.append(icon("file", bx + 52, y + 119, "acc", 24))
    parts.append(text(bx + 92, y + 112, "05_niche.h5ad", 18, 700, YELLOW_DARK))
    parts.append(text(bx + 92, y + 136, u["niche_out"], 13, 400, YELLOW_DARK))
    y += 252
    parts.append(arrow_down(W / 2, y + 2, y + 30, PALETTE["acc"][0]))
    y += 34

    # 下游三张卡片
    efg = p["efg"]
    rows = max(len(efg["E"]), len(efg["F"]) or 1, len(efg["G"]))
    card_h = 92 + rows * 78 + (44 if efg["na"] else 0)
    parts.append(rect(M - 8, y, W - 2 * M + 16, card_h + 58, fill=YELLOW[1], stroke=PALETTE["acc"][0], rx=18, sw=1.5, dash="7 5"))
    parts.append(text(M + 14, y + 30, u["downstream"], 14, 600, YELLOW_DARK))
    cw = (W - 2 * M - 40) / 3
    for j, key in enumerate("EFG"):
        cx0, cy0 = M + 10 + j * (cw + 10), y + 44
        parts.append(rect(cx0, cy0, cw, card_h, shadow=True))
        parts.append(f'<path d="M{cx0 + 16},{cy0} h{cw - 32} a16,16 0 0 1 16,16 v14 h{-cw} v-14 a16,16 0 0 1 16,-16 z" fill="{PALETTE[key][1]}"/>')
        parts.append(icon({"E": "nhood", "F": "comm", "G": "special"}[key], cx0 + 38, cy0 + 44, key, 22))
        parts.append(text(cx0 + 70, cy0 + 50, u[key], 16, 700, PALETTE[key][0]))
        items = efg[key]
        if not items:
            parts.append(text(cx0 + 24, cy0 + 100, u["na"], 14, 600, SUB))
            parts.append(text(cx0 + 24, cy0 + 120, "no ligand–receptor pairs in antibody panels" if lang == "en" else "抗体 panel 里基本没有配受体", 12, 400, SUB))
        for i, st in enumerate(items):
            iy = cy0 + 82 + i * 78
            code = f"{key}{i + 1}"
            svg, w = pill(cx0 + 20, iy, code, key, 12, solid=not st["optional"], dashed=st["optional"])
            parts.append(svg)
            parts.append(text(cx0 + 30 + w, iy + 16, st[lang], 15, 600, TEXT))
            sub = st["script"] + (f" · {u['optional']}" if st["optional"] else "")
            parts.append(text(cx0 + 30 + w, iy + 36, sub, 12, 400, SUB))
            uses_en, uses_zh, out = SCRIPT_IO[st["script"]]
            parts.append(text(cx0 + 30 + w, iy + 54, ("读 " + uses_zh) if lang == "zh" else ("reads " + uses_en), 11, 400, SUB))
            parts.append(text(cx0 + 30 + w, iy + 69, ("写 " if lang == "zh" else "writes ") + out, 11, 600, YELLOW_DARK))
        if key == "G" and efg["na"]:
            parts.append(text(cx0 + 20, cy0 + card_h - 34, f"{u['na']}:", 12, 600, SUB))
            parts.append(text(cx0 + 20, cy0 + card_h - 16, efg["na"][0 if lang == "en" else 1], 12, 400, SUB))
    y += card_h + 58 + 22

    # 溯源
    parts.append(rect(M - 8, y, W - 2 * M + 16, 76, fill=GREY[1], stroke=BORDER, rx=16, dash="6 5"))
    parts.append(icon("prov", M + 34, y + 38, "Z", 22))
    parts.append(text(M + 70, y + 32, u["prov"], 15, 700, PALETTE["Z"][0]))
    parts.append(text(M + 70, y + 54, u["prov_sub"], 12, 400, SUB))
    y += 96
    return wrap(W, y, "".join(parts))

def wrap(w, h, body):
    defs = ('<defs><filter id="shadow" x="-5%" y="-5%" width="110%" height="120%">'
            f'<feDropShadow dx="0" dy="2" stdDeviation="3" flood-color="{GREY_DARK}" flood-opacity="0.07"/></filter></defs>')
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h:.0f}" viewBox="0 0 {w} {h:.0f}" font-family="{FONT}">'
            f'{defs}<rect width="100%" height="100%" fill="{CANVAS}"/>{body}</svg>\n')

# %%
# 步骤：通用小部件（折线箭头、菱形判断框、自动换行）
# 输入：坐标、文字
# 输出：SVG 片段
def poly_arrow(points, color=GREY[0], width=2):
    d = "M" + " L".join(f"{x:.1f},{y:.1f}" for x, y in points)
    (x1, y1), (x2, y2) = points[-2], points[-1]
    if x1 == x2:
        sign = 1 if y2 > y1 else -1
        head = f"M{x2 - 6},{y2 - 10 * sign} L{x2},{y2} L{x2 + 6},{y2 - 10 * sign} Z"
    else:
        sign = 1 if x2 > x1 else -1
        head = f"M{x2 - 10 * sign},{y2 - 6} L{x2},{y2} L{x2 - 10 * sign},{y2 + 6} Z"
    return f'<path d="{d}" fill="none" stroke="{color}" stroke-width="{width}" stroke-linejoin="round"/><path d="{head}" fill="{color}"/>'

def diamond(cx, cy, w, h, key, lines):
    color, light = PALETTE[key]
    s = f'<path d="M{cx},{cy - h / 2} L{cx + w / 2},{cy} L{cx},{cy + h / 2} L{cx - w / 2},{cy} Z" fill="{light}" stroke="{color}" stroke-width="1.8"/>'
    y0 = cy - (len(lines) - 1) * 9 + 5
    return s + "".join(text(cx, y0 + i * 18, line, 14, 600, TEXT, "middle") for i, line in enumerate(lines))

def wrap_text(s, size, max_width):
    if width_of(s, size) <= max_width:
        return [s]
    if " · " in s or " " in s:
        sep = " · " if " · " in s else " "
        words, lines, line = s.split(sep), [], ""
        for w in words:
            trial = (line + sep + w) if line else w
            if width_of(trial, size) > max_width and line:
                lines.append(line); line = w
            else:
                line = trial
        return lines + [line]
    lines, line = [], ""
    for c in s:
        if width_of(line + c, size) > max_width:
            lines.append(line); line = c
        else:
            line += c
    return lines + [line]

LANES = {   # 数据类型：图标、模块色、英文名、中文名、平台、对应流程
    "spot": ("data", "A", "Spot-level", "spot 级", "Visium · Stereo-seq bin≥50", ("P1", "P2")),
    "hd": ("domain", "B", "High-res sequencing", "近单细胞，测序型", "Visium HD · Stereo-seq bin≤20 · Slide-seq", ("P3", "P4")),
    "imaging": ("cluster", "C", "Imaging, targeted RNA", "单细胞，靶向 RNA", "Xenium · CosMx · MERSCOPE", ("P5", "P6")),
    "protein": ("nhood", "G", "Spatial proteomics", "空间蛋白组", "CODEX · IMC · mIHC · MIBI", ("P7",)),
}

LANE_COLOR = {"P1": "A", "P2": "A", "P3": "B", "P4": "B", "P5": "C", "P6": "C", "P7": "G"}   # 流程编号按数据类型着色

def ref_label(pid, zh):
    if pid == "P7":
        return "参考不适用" if zh else "reference n/a"
    return ("有单细胞参考" if zh else "with reference") if int(pid[1]) % 2 == 1 else ("无单细胞参考" if zh else "no reference")

# %%
# 步骤：画决策流程图（从数据出发，一步步问到 P1–P7）
# 输入：语言
# 输出：SVG 字符串
def draw_decision(lang):
    zh = lang == "zh"
    W = 1460
    by_id = {p["id"]: p for p in PIPELINES}
    parts = [text(36, 52, "图 1 · 选择流程：我的数据用哪一套？" if zh else "Figure 1 · Choose a pipeline: which one fits my data?", 26, 700),
             text(36, 80, "按顺序回答三个问题，落到 P1–P7 中的一套 → 图 2 看这套流程每一步用什么 → 图 3 打开这套流程的详细流程图" if zh else
                  "Answer three questions to land on one of P1–P7 → Figure 2 shows what it runs at each step → Figure 3 is the detailed diagram of that pipeline", 14, 400, SUB)]
    # 起点
    parts.append(rect(600 - 150, 104, 300, 58, fill=PALETTE["input"][1], stroke=PALETTE["input"][0], rx=29))
    parts.append(text(600, 131, "你的空间组学数据" if zh else "Your spatial data", 16, 700, TEXT, "middle"))
    parts.append(text(600, 151, "从仪器原始输出开始" if zh else "start from the instrument output", 12, 400, SUB, "middle"))
    parts.append(poly_arrow([(600, 162), (600, 190)]))
    # 问题 1
    parts.append(diamond(600, 244, 320, 104, "Z", ["问题 1：一个观测点", "是不是一个细胞？"] if zh else ["Q1 · is one observation", "one cell?"]))
    parts.append(poly_arrow([(440, 244), (225, 244), (225, 480)]))
    parts.append(text(330, 232, "否：spot 里混有多个细胞" if zh else "No · each spot mixes cells", 13, 600, PALETTE["A"][0], "middle"))
    parts.append(poly_arrow([(760, 244), (960, 244), (960, 316)]))
    parts.append(text(860, 232, "是：单细胞级" if zh else "Yes · single-cell level", 13, 600, PALETTE["B"][0], "middle"))
    # 问题 2
    parts.append(diamond(960, 362, 280, 92, "Z", ["问题 2：测的是什么？"] if zh else ["Q2 · what is measured?"]))
    parts.append(f'<path d="M960,408 L960,436 M625,436 L1300,436" fill="none" stroke="{GREY[0]}" stroke-width="2"/>')
    for cx, label, key in [(625, "全转录组（测序）" if zh else "whole transcriptome (sequencing)", "B"),
                           (1015, "靶向 RNA panel（成像）" if zh else "targeted RNA panel (imaging)", "C"),
                           (1300, "蛋白 panel" if zh else "protein panel", "G")]:
        parts.append(poly_arrow([(cx, 436), (cx, 480)]))
        parts.append(text(cx + 8, 462, label, 12, 600, PALETTE[key][0]))
    # 数据类型框
    for lane, cx in [("spot", 225), ("hd", 625), ("imaging", 1015), ("protein", 1300)]:
        kind, key, en, zh_name, platforms, pids = LANES[lane]
        parts.append(rect(cx - 135, 480, 270, 82, fill=PALETTE[key][1], stroke=PALETTE[key][0], rx=14, sw=1.5))
        parts.append(icon(kind, cx - 100, 521, key, 20))
        parts.append(text(cx - 70, 515, zh_name if zh else en, 15, 700, TEXT))
        for i, line in enumerate(wrap_text(platforms, 11, 196)):
            parts.append(text(cx - 70, 535 + i * 15, line, 11, 400, SUB))
    # 问题 3 和叶子
    leaf_x = {"P1": 130, "P2": 320, "P3": 530, "P4": 720, "P5": 920, "P6": 1110, "P7": 1300}
    for cx, left, right in [(225, "P1", "P2"), (625, "P3", "P4"), (1015, "P5", "P6")]:
        parts.append(poly_arrow([(cx, 562), (cx, 604)]))
        parts.append(diamond(cx, 650, 250, 92, "C", ["问题 3：有带原始 counts", "的单细胞参考吗？"] if zh else ["Q3 · scRNA-seq reference", "with raw counts?"]))
        parts.append(f'<path d="M{cx},696 L{cx},716 M{leaf_x[left]},716 L{leaf_x[right]},716" fill="none" stroke="{GREY[0]}" stroke-width="2"/>')
        for pid, label in [(left, "有" if zh else "yes"), (right, "无" if zh else "no")]:
            parts.append(poly_arrow([(leaf_x[pid], 716), (leaf_x[pid], 760)]))
            parts.append(text(leaf_x[pid] + 8, 743, label, 12, 600, PALETTE["C"][0]))
    parts.append(poly_arrow([(1300, 562), (1300, 760)]))
    parts.append(text(1290, 690, "单细胞参考不适用" if zh else "reference not applicable", 12, 600, PALETTE["G"][0], "end"))
    for pid, cx in leaf_x.items():
        p = by_id[pid]
        key = LANE_COLOR[pid]
        with_ref = pid in ("P1", "P3", "P5")
        parts.append(rect(cx - 92, 760, 184, 112, shadow=True))
        svg, pw = pill(cx - 78, 774, pid, key, 14, solid=True)
        parts.append(svg)
        svg, _ = pill(cx - 78 + pw + 8, 777, ref_label(pid, zh), "acc" if with_ref else "input", 10)
        parts.append(svg)
        parts.append(text(cx - 78, 826, p["annot"][lang], 14, 700, TEXT))
        if p["annot"]["lang"] == "R":     # 只在用到 R 的方法旁标一个小 R
            svg, _ = pill(cx - 78 + width_of(p["annot"][lang], 14) + 6, 812, "R", "R", 9, solid=True)
            parts.append(svg)
        file = {"P1": "P1_lowres_scrna.sh", "P2": "P2_lowres_noref.sh", "P3": "P3_hd_scrna.sh", "P4": "P4_hd_noref.sh",
                "P5": "P5_imaging_scrna.sh", "P6": "P6_imaging_noref.sh", "P7": "P7_protein.sh"}[pid]
        parts.append(text(cx - 78, 848, file, 11, 400, SUB))
        parts.append(text(cx - 78, 863, p["annot"]["script"], 10.5, 400, GREY[0]))
    return wrap(W, 896, "".join(parts))

# %%
# 步骤：分步对照图的内容：每一行是一个步骤，每一列是一种数据类型
# 格子：(状态, 英文标题, 中文标题, 英文说明, 中文说明, 是否 R 步骤)；状态 run = 默认运行，opt = 可选，na = 不适用
def cell(state, en, zh, sub_en="", sub_zh="", r=False):
    return dict(state=state, en=en, zh=zh, sub_en=sub_en, sub_zh=sub_zh, r=r)

ROW_IO = {'Check config': ('reads config.yaml, samples.tsv', '读 config.yaml、samples.tsv', 'writes nothing (check only)', '写 无（只检查）'), 'Load': ('reads instrument output', '读 仪器原始输出', 'writes 01_raw.h5ad', '写 01_raw.h5ad'), 'QC & clustering': ('reads 01_raw · layers counts', '读 01_raw · layers counts', 'writes 02_clustered.h5ad', '写 02_clustered.h5ad'), 'Annotation|with': ('reads 02_clustered · counts + ref', '读 02_clustered · counts + 参考', 'writes 03_annotated.h5ad', '写 03_annotated.h5ad'), 'Annotation|without': ('reads 02_clustered · X + markers', '读 02_clustered · X + markers', 'writes 03_annotated.h5ad', '写 03_annotated.h5ad'), 'Spatial domain': ('reads 03_annotated · X_pca + coords', '读 03_annotated · X_pca + 坐标', 'writes 04_domain.h5ad', '写 04_domain.h5ad'), 'Niche': ('reads 04_domain · types or composition', '读 04_domain · 细胞类型或组成', 'writes 05_niche.h5ad', '写 05_niche.h5ad'), 'Neighbourhoods': ('reads 05_niche · cell_type / composition', '读 05_niche · cell_type 或组成', 'writes nhood / misty tables', '写 邻域、MISTy 图表'), 'Communication': ('reads 05_niche · X + coords', '读 05_niche · X + 坐标', 'writes lr_*.csv', '写 lr_*.csv'), 'Spatial CNV': ('reads 05_niche · X + cell_type', '读 05_niche · X + cell_type', 'writes cnv_per_spot.csv', '写 cnv_per_spot.csv'), 'Spatially variable': ('reads 05_niche · X + coords', '读 05_niche · X + 坐标', 'writes svg_moran.csv', '写 svg_moran.csv'), 'NMF programmes': ('reads 05_niche · X (HVGs)', '读 05_niche · X（高变基因）', 'writes nmf_*.csv', '写 nmf_*.csv'), 'Niche subtypes': ('reads 05_niche · nhood_composition', '读 05_niche · nhood_composition', 'writes subtype_*.csv', '写 subtype_*.csv'), 'Hotspots & trends': ('reads 05_niche · composition + coords', '读 05_niche · 组成 + 坐标', 'writes hotspot / tendency csv', '写 hotspot、tendency 表'), 'Provenance': ('reads config + all results', '读 config 和全部结果', 'writes results/provenance/', '写 results/provenance/')}

# 每一行的样本处理方式：(类型, 英文标签, 中文标签)；类型 joint 联合、per 逐样本、single 单样本、needs 需要多样本
ROW_MODE = {
    "Check config": [],
    "Load": [("per", "per sample → merged", "逐样本读入 → 合并")],
    "QC & clustering": [("joint", "joint (Harmony)", "多样本联合（Harmony）")],
    "Annotation|with": [("per", "RCTD per sample", "RCTD 逐样本"), ("joint", "transfer joint", "标签转移联合")],
    "Annotation|without": [("joint", "joint", "多样本联合")],
    "Spatial domain": [("joint", "joint", "多样本联合")],
    "Niche": [("joint", "joint", "多样本联合")],
    "Neighbourhoods": [("per", "per sample → mean", "逐样本 → 平均"), ("needs", "diff: ≥2/group", "差异需多样本")],
    "Communication": [("per", "LIANA+ per sample", "LIANA+ 逐样本"), ("single", "COMMOT single", "COMMOT 单样本")],
    "Spatial CNV": [("joint", "joint", "多样本联合")],
    "Spatially variable": [("per", "per sample → summary", "逐样本 → 汇总")],
    "NMF programmes": [("joint", "joint", "多样本联合")],
    "Niche subtypes": [("joint", "joint", "多样本联合")],
    "Hotspots & trends": [("per", "hotspot per sample", "热点逐样本"), ("single", "trend single", "趋势单样本")],
    "Provenance": [],
}
MODE_STYLE = {"joint": ("Py", False), "per": ("input", False), "single": ("input", True), "needs": ("acc", False)}   # (配色, 是否虚线)

STEP_ROWS = [
    dict(key="Z", en="Check config", zh="检查配置", sub_en="validate_config.py", sub_zh="validate_config.py",
         span=cell("run", "validate_config.py", "validate_config.py", "checks config.yaml against the schema and that every sample path exists", "按 schema 检查 config.yaml，并确认每个样本路径都存在")),
    dict(key="A", en="Load", zh="读入", sub_en="→ 01_raw.h5ad, coordinates in µm", sub_zh="→ 01_raw.h5ad，坐标换算成 µm", cells=[
        cell("run", "load_visium.py", "load_visium.py", "Stereo-seq: load_stereo.py", "Stereo-seq 用 load_stereo.py"),
        cell("run", "load_visium_hd.py", "load_visium_hd.py", "16 µm bins; cellbin: load_stereo.py", "16 µm bin；cellbin 用 load_stereo.py"),
        cell("run", "load_xenium.py", "load_xenium.py", "CosMx: load_cosmx.py (global px)", "CosMx 用 load_cosmx.py（全片坐标）"),
        cell("run", "load_codex.py", "load_codex.py", "cell × marker intensity CSV", "细胞 × marker 强度表")]),
    dict(key="B", en="QC & clustering", zh="质控与聚类", sub_en="→ 02_clustered.h5ad", sub_zh="→ 02_clustered.h5ad", cells=[
        cell("run", "preprocess_rna.py", "preprocess_rna.py", "min_counts 500 · 2000 HVGs · Harmony", "min_counts 500 · 2000 高变基因 · Harmony"),
        cell("run", "preprocess_rna.py", "preprocess_rna.py", "min_counts 50 · 2000 HVGs · Harmony", "min_counts 50 · 2000 高变基因 · Harmony"),
        cell("run", "preprocess_rna.py", "preprocess_rna.py", "all panel genes, no HVG selection", "panel 基因全部使用，不选高变基因"),
        cell("run", "preprocess_protein.py", "preprocess_protein.py", "arcsinh → z-score → Leiden", "arcsinh → z-score → Leiden")]),
    dict(key="C", en="Annotation", zh="注释", sub_en="with ref", sub_zh="有单细胞参考时", cells=[
        cell("run", "RCTD full · P1", "RCTD full · P1", "deconvolution → cell_props; alt. cell2location", "去卷积 → cell_props；也可用 cell2location", True),
        cell("run", "RCTD doublet · P3", "RCTD doublet · P3", "segmented cells: label transfer", "细胞分割结果可改用标签转移", True),
        cell("run", "Label transfer · P5", "标签转移 · P5", "Harmony + KNN → cell_type_score", "Harmony + KNN → cell_type_score"),
        cell("na", "Not applicable", "不适用", "scRNA reference does not fit protein", "RNA 参考不适用于蛋白数据")]),
    dict(key="C", en="Annotation", zh="注释", sub_en="without reference", sub_zh="没有单细胞参考时", cells=[
        cell("run", "Marker scores · P2", "marker 打分 · P2", "cell_scores (not proportions)", "cell_scores（不是比例）"),
        cell("run", "Marker annotation · P4", "marker 注释 · P4", "per Leiden cluster", "按 Leiden 聚类整簇注释"),
        cell("run", "Marker annotation · P6", "marker 注释 · P6", "panel genes only", "markers.yaml 只写 panel 里的基因"),
        cell("run", "Marker annotation · P7", "marker 注释 · P7", "per cluster + manual check", "按聚类注释，再人工核对")]),
    dict(key="D", en="Spatial domain", zh="空间域", sub_en="→ 04_domain.h5ad", sub_zh="→ 04_domain.h5ad",
         span=cell("run", "domain_cellcharter.py", "domain_cellcharter.py", "same for all data types · or BANKSY in R (domain_banksy.Rmd, λ = 0.8)", "所有数据类型相同 · 也可用 R 的 BANKSY（domain_banksy.Rmd，λ = 0.8）")),
    dict(key="D", en="Niche", zh="niche", sub_en="→ 05_niche.h5ad", sub_zh="→ 05_niche.h5ad", cells=[
        cell("run", "Cell niche", "细胞 niche", "cluster cell proportions of each spot", "对每个 spot 的细胞比例聚类"),
        cell("run", "Niche / RCN", "niche / RCN", "neighbour composition within radius", "统计半径内邻居的细胞组成"),
        cell("run", "RCN", "RCN", "neighbour composition within radius", "统计半径内邻居的细胞组成"),
        cell("run", "CN", "CN", "cellular neighbourhoods", "细胞邻域（CN）")]),
    dict(key="E", en="Neighbourhoods", zh="邻域与共定位", sub_en="co-localisation, group differences", sub_zh="共定位、组间差异", cells=[
        cell("run", "MISTy co-localisation", "MISTy 共定位", "no nhood enrichment: spots mix cells", "不用邻域富集：spot 里混有多个细胞"),
        cell("run", "Nhood enrichment + group diff", "邻域富集 + 组间差异", "diff needs ≥ 2 groups × ≥ 2 samples", "组间差异要求 ≥ 2 组、每组 ≥ 2 个样本"),
        cell("run", "Nhood enrichment + group diff", "邻域富集 + 组间差异", "diff needs ≥ 2 groups × ≥ 2 samples", "组间差异要求 ≥ 2 组、每组 ≥ 2 个样本"),
        cell("run", "Nhood enrichment, diff, MISTy", "邻域富集、组间差异、MISTy", "all three run by default", "三项默认都运行")]),
    dict(key="F", en="Communication", zh="细胞通讯", sub_en="ligand–receptor", sub_zh="配体-受体", cells=[
        cell("run", "LIANA+ spatial LR", "LIANA+ 空间配受体", "optional: COMMOT signal flow", "可选：COMMOT 信号流向"),
        cell("run", "LIANA+ spatial LR", "LIANA+ 空间配受体", "optional: LR within neighbourhoods", "可选：邻域内配受体"),
        cell("run", "LIANA+ spatial LR", "LIANA+ 空间配受体", "check panel LR coverage first", "先看 panel 覆盖了多少配受体对"),
        cell("na", "Not applicable", "不适用", "no ligand–receptor pairs in antibody panels", "抗体 panel 里基本没有配受体")]),
    dict(key="G", en="Spatial CNV", zh="空间 CNV", sub_en="tumour samples", sub_zh="肿瘤样本", cells=[
        cell("opt", "cnv_infercnvpy.py", "cnv_infercnvpy.py", "needs GTF + normal reference types", "需要 GTF 和正常参考细胞类型"),
        cell("opt", "cnv_infercnvpy.py", "cnv_infercnvpy.py", "use 16 µm bins or segmented cells", "用 16 µm bin 或细胞分割结果"),
        cell("na", "Not applicable", "不适用", "panel too small for CNV", "panel 基因太少"),
        cell("na", "Not applicable", "不适用", "no genome positions", "蛋白没有基因组位置")]),
    dict(key="G", en="Spatially variable", zh="空间可变特征", sub_en="Moran's I", sub_zh="Moran's I",
         span=cell("run", "svg_moran.py", "svg_moran.py", "genes for RNA data, markers for protein data", "RNA 数据找空间可变基因，蛋白数据找空间可变 marker")),
    dict(key="G", en="NMF programmes", zh="NMF 表达程序", sub_en="gene programmes", sub_zh="基因程序", cells=[
        cell("run", "nmf_programs.py", "nmf_programs.py", "programmes per niche", "每个 niche 的程序使用量"),
        cell("run", "nmf_programs.py", "nmf_programs.py", "fit on a subsample", "先下采样再拟合"),
        cell("na", "Not in pipeline", "不在流程中", "too few genes for programmes", "基因太少，拆不出程序"),
        cell("na", "Not in pipeline", "不在流程中", "too few markers", "marker 太少")]),
    dict(key="G", en="Niche subtypes", zh="微环境亚型", sub_en="NMF on neighbourhoods", sub_zh="对邻域组成做 NMF", cells=[
        cell("na", "Not applicable", "不适用", "spots are not single cells", "spot 不是单个细胞"),
        cell("opt", "nmf_subtype.py", "nmf_subtype.py", "set target_cell_type first", "先设定 target_cell_type"),
        cell("opt", "nmf_subtype.py", "nmf_subtype.py", "set target_cell_type first", "先设定 target_cell_type"),
        cell("opt", "nmf_subtype.py", "nmf_subtype.py", "set target_cell_type first", "先设定 target_cell_type")]),
    dict(key="G", en="Hotspots & trends", zh="热点与空间趋势", sub_en="optional for all", sub_zh="都可选",
         span=cell("opt", "hotspot_getis.py · spatial_tendency_soapy.py", "hotspot_getis.py · spatial_tendency_soapy.py", "Getis-Ord hotspots and distances; gene / cell trends away from a region", "Getis-Ord 热点和热点间距离；离某个区域越远，基因或细胞如何变化")),
    dict(key="Z", en="Provenance", zh="溯源记录", sub_en="last step of every pipeline", sub_zh="每套流程的最后一步",
         span=cell("run", "record_provenance.py", "record_provenance.py", "results/provenance/: config copy, software versions, SHA-256 of inputs and outputs", "results/provenance/：配置副本、软件版本、输入和输出文件的 SHA-256")),
]

# %%
# 步骤：画分步对照图（泳道：每种数据类型一列，从上到下是步骤）
# 输入：语言
# 输出：SVG 字符串
def draw_cell(x, y, w, h, c, key, lang):
    color, light = PALETTE[key]
    if c["state"] == "run":
        s = rect(x, y, w, h, fill=light, stroke=color, rx=12, sw=1.2)
        title_color = TEXT
    elif c["state"] == "opt":
        s = rect(x, y, w, h, fill=CANVAS, stroke=color, rx=12, sw=1.2, dash="5 4")
        title_color = TEXT
    else:
        s = rect(x, y, w, h, fill=GREY[1], stroke=GREY_LINE, rx=12, sw=1.2)
        title_color = GREY[0]
    title = ("✕ " if c["state"] == "na" else "") + c[lang]
    s += text(x + 14, y + 24, title, 13.5, 700, title_color)
    tag_x = x + w - 12
    if c["r"]:
        svg, pw = pill(tag_x - (width_of("R", 10) + 18), y + 9, "R", "R", 10, solid=True)
        s += svg; tag_x -= pw + 6
    if c["state"] == "opt":
        label = "可选" if lang == "zh" else "optional"
        svg, _ = pill(tag_x - (width_of(label, 10) + 18), y + 9, label, key, 10, dashed=True)
        s += svg
    for i, line in enumerate(wrap_text(c["sub_" + lang], 11.5, w - 40)[:2]):
        s += text(x + 14, y + 44 + i * 15, line, 11.5, 400, SUB if c["state"] != "na" else GREY[0])
    return s

def draw_steps(lang):
    zh = lang == "zh"
    W, M, label_w, gap = 1440, 30, 250, 8
    lane_w = (W - 2 * M - label_w - 3 * gap) / 4
    lane_x = [M + label_w + i * (lane_w + gap) for i in range(4)]
    parts = [text(M + 6, 52, "图 2 · 分步总览：七套流程每一步用什么" if zh else "Figure 2 · Step-by-step overview of all seven pipelines", 26, 700),
             text(M + 6, 80, "每一列是一种数据（表头标出对应的流程），从上往下是运行顺序；左侧标出每一步读写的文件和样本处理方式；单套流程的细节见图 3" if zh else
                  "Each column is one data type (pipelines named in the header); rows run top to bottom; the left column gives each step's files and how it handles samples. Single-pipeline detail: Figure 3", 14, 400, SUB)]
    x = M + 6
    for label, kind in [("默认运行" if zh else "runs by default", "run"), ("可选" if zh else "optional", "opt"), ("不适用" if zh else "not applicable", "na")]:
        w = width_of(label, 12) + 34
        fill, stroke, dash = {"run": (PALETTE["C"][1], PALETTE["C"][0], None), "opt": (CANVAS, PALETTE["C"][0], "4 3"), "na": (GREY[1], GREY_LINE, None)}[kind]
        parts.append(rect(x, 96, w, 24, fill=fill, stroke=stroke, rx=8, sw=1.2, dash=dash))
        parts.append(text(x + w / 2, 113, label, 12, 600, TEXT if kind != "na" else GREY[0], "middle"))
        x += w + 10
    svg, w = pill(x, 97, "R", "R", 11, solid=True)
    parts.append(svg + text(x + w + 6, 113, "R 步骤（其余都是 Python）" if zh else "R step (all others are Python)", 12, 400, SUB))
    x += w + 6 + width_of("R 步骤（其余都是 Python）" if zh else "R step (all others are Python)", 12) + 30
    parts.append(text(x, 113, "样本处理：" if zh else "Samples:", 12, 600, SUB))
    x += width_of("样本处理：" if zh else "Samples:", 12) + 8
    for kind, label in [("joint", "多样本联合" if zh else "joint"), ("per", "逐样本" if zh else "per sample"),
                        ("single", "单样本" if zh else "single sample"), ("needs", "需多样本" if zh else "needs ≥2 samples")]:
        color_key, dashed = MODE_STYLE[kind]
        svg, tw = pill(x, 99, label, color_key, 10, dashed=dashed)
        parts.append(svg); x += tw + 6
    y = 140
    # 泳道表头
    for i, lane in enumerate(LANES):
        kind, key, en, zh_name, platforms, pids = LANES[lane]
        parts.append(rect(lane_x[i], y, lane_w, 104, fill=PALETTE[key][1], stroke=PALETTE[key][0], rx=14, sw=1.5))
        parts.append(icon(kind, lane_x[i] + 36, y + 40, key, 20))
        parts.append(text(lane_x[i] + 66, y + 34, zh_name if zh else en, 15, 700, TEXT))
        for j, line in enumerate(wrap_text(platforms, 11, lane_w - 80)):
            parts.append(text(lane_x[i] + 66, y + 53 + j * 14, line, 11, 400, SUB))
        px = lane_x[i] + 16
        for pid, ref in zip(pids, ([("有参考" if zh else "with ref"), ("无参考" if zh else "no ref")] if len(pids) == 2 else [("参考不适用" if zh else "no ref needed")])):
            svg, pw = pill(px, y + 74, f"{pid} · {ref}", key, 11, solid=pid in ("P1", "P3", "P5", "P7"))
            parts.append(svg); px += pw + 8
    y += 104 + 22
    # 步骤行
    row_h, step_no, last_title = 96, 0, None
    for row in STEP_ROWS:
        key = row["key"]
        if row[lang] != last_title:
            step_no += 1
        new_step = row[lang] != last_title
        last_title = row[lang]
        color, light = PALETTE[key]
        if new_step:
            parts.append(f'<circle cx="{M + 18}" cy="{y + 20}" r="14" fill="{color}"/>' + text(M + 18, y + 25, str(step_no), 13, 700, "#FFFFFF", "middle"))
            parts.append(text(M + 42, y + 25, row[lang], 15, 700, TEXT))
            if row["key"] == "C":                    # 注释分两行：有参考 / 无参考
                parts.append(text(M + 50 + width_of(row[lang], 15), y + 25, row["sub_" + lang], 12, 600, SUB))
        else:
            parts.append(text(M + 42, y + 25, row["sub_" + lang], 12, 600, SUB))
        io_key = row["en"] + ("|without" if row["sub_en"].startswith("without") else "|with" if row["sub_en"] == "with ref" else "")
        reads_en, reads_zh, writes_en, writes_zh = ROW_IO[io_key]
        parts.append(text(M + 42, y + 43, reads_zh if zh else reads_en, 10.5, 400, SUB))
        parts.append(text(M + 42, y + 58, writes_zh if zh else writes_en, 10.5, 600, YELLOW_DARK))
        tag_x, tag_y = M + 42, y + 66
        for kind, en, zh_label in ROW_MODE[io_key]:
            color_key, dashed = MODE_STYLE[kind]
            label = zh_label if zh else en
            if tag_x + width_of(label, 9) + 18 > M + label_w - 4:        # 放不下就换到下一行
                tag_x, tag_y = M + 42, tag_y + 22
            svg, tw = pill(tag_x, tag_y, label, color_key, 9, dashed=dashed)
            parts.append(svg); tag_x += tw + 4
        if "span" in row:
            span_w = lane_x[3] + lane_w - lane_x[0]
            parts.append(draw_cell(lane_x[0], y, span_w, row_h - 12, row["span"], key, lang))
        else:
            for i, c in enumerate(row["cells"]):
                parts.append(draw_cell(lane_x[i], y, lane_w, row_h - 12, c, key, lang))
        if row is not STEP_ROWS[-1]:
            for i in range(4):
                cx = lane_x[i] + lane_w / 2
                parts.append(poly_arrow([(cx, y + row_h - 11), (cx, y + row_h)], GREY[0], 1.6))
        y += row_h
    return wrap(W, y + 20, "".join(parts))

# %%
# 步骤：图 4 的内容：5 个阶段文件，每个文件由谁写出、读上一步的什么、新增哪些字段
# 输入：无
# 输出：STAGES、CONSUMERS
STAGES = [
    dict(file="01_raw.h5ad", by="load_*.py", reads_en="instrument output", reads_zh="仪器原始输出",
         adds=[("layers['counts']", "raw counts", "原始 counts（蛋白为 layers['raw']）"), ("obsm['spatial']", "native coordinates", "平台原始坐标"),
               ("obsm['spatial_um']", "coordinates in µm", "µm 坐标"), ("obs sample / condition", "sample and group", "样本和分组")]),
    dict(file="02_clustered.h5ad", by="preprocess_*.py", reads_en="01_raw · layers['counts']", reads_zh="01_raw · layers['counts']",
         adds=[("X", "log-normalised (protein: z-score)", "log 标准化（蛋白：z-score）"), ("var highly_variable", "genes used for PCA", "PCA 用的基因"),
               ("obsm X_pca(_harmony)", "embedding", "降维结果"), ("obs['leiden']", "expression clusters", "表达聚类")]),
    dict(file="03_annotated.h5ad", by="C_annotate/*", reads_en="02_clustered · counts (RCTD, transfer) or X (markers)", reads_zh="02_clustered · counts（RCTD、标签转移）或 X（marker）",
         adds=[("obs['cell_type']", "cell type / dominant type", "细胞类型 / 主导类型"), ("obsm['cell_props']", "proportions (RCTD, c2l)", "比例（RCTD、cell2location）"),
               ("obsm['cell_scores']", "marker scores", "marker 打分"), ("obs cell_type_score", "label-transfer confidence", "标签转移置信度")]),
    dict(file="04_domain.h5ad", by="domain_cellcharter.py", reads_en="03_annotated · X_pca + spatial_um", reads_zh="03_annotated · X_pca + spatial_um",
         adds=[("obs['domain']", "spatial domain", "空间域"), ("obsm X_cellcharter", "neighbour-aggregated embedding", "聚合邻居后的降维")]),
    dict(file="05_niche.h5ad", by="niche_composition.py", reads_en="04_domain · cell_type or cell_props + spatial_um", reads_zh="04_domain · cell_type 或 cell_props + spatial_um",
         adds=[("obs['niche']", "niche / RCN / CN", "niche / RCN / CN"), ("obsm nhood_composition", "neighbourhood composition", "邻域细胞组成")]),
]
CONSUMERS = [
    ("X + obsm['spatial_um']", "log expression and coordinates", "log 表达和 µm 坐标",
     ["lr_spatial_liana", "svg_moran", "nmf_programs", "commot_flow", "cnv_infercnvpy", "spatial_tendency_soapy"]),
    ("obs['cell_type'] + spatial_um", "cell types and coordinates", "细胞类型和坐标",
     ["nhood_enrichment", "diff_nhood", "ligrec_neighborhood"]),
    ("obsm cell_props / cell_scores / nhood_composition", "cell-type composition", "细胞组成",
     ["misty_colocalization", "hotspot_getis", "nmf_subtype"]),
    ("obs['niche'] / obs['domain']", "grouping for summaries and plots", "按 niche / 空间域分组汇总、画图",
     ["nmf_programs", "spatial_tendency_soapy", "all figures"]),
]

# %%
# 步骤：画图 4（数据流）
# 输入：语言
# 输出：SVG 字符串
def draw_dataflow(lang):
    zh = lang == "zh"
    W, M, card_w, gap = 1440, 30, 252, 25
    parts = [text(M + 6, 52, "图 4 · 数据流：每一步读什么、写什么" if zh else "Figure 4 · Data flow: what each step reads and writes", 26, 700),
             text(M + 6, 80, "七套流程共用同一条数据链：每一步读上一步文件里的某一部分，在它的基础上加字段，写出下一个文件；下游脚本都读 05_niche.h5ad" if zh else
                  "All seven pipelines share one data chain: each step reads part of the previous file, adds fields and writes the next file; downstream scripts all read 05_niche.h5ad", 14, 400, SUB)]
    y = 108
    for i, st in enumerate(STAGES):
        x = M + i * (card_w + gap)
        parts.append(rect(x, y, card_w, 262, shadow=True))
        parts.append(rect(x, y, card_w, 54, fill=YELLOW[1], stroke=YELLOW[0], rx=16, sw=1.5))
        parts.append(icon("file", x + 28, y + 27, "acc", 16))
        parts.append(text(x + 52, y + 25, st["file"], 15, 700, YELLOW_DARK))
        parts.append(text(x + 52, y + 43, ("写出：" if zh else "written by ") + st["by"], 11, 400, YELLOW_DARK))
        parts.append(text(x + 14, y + 76, "读上一步" if zh else "reads from previous step", 11, 600, SUB))
        for j, line in enumerate(wrap_text(st["reads_" + lang], 11, card_w - 28)[:2]):
            parts.append(text(x + 14, y + 92 + j * 14, line, 11, 400, TEXT))
        parts.append(f'<line x1="{x + 14}" y1="{y + 118}" x2="{x + card_w - 14}" y2="{y + 118}" stroke="{BORDER}" stroke-width="1"/>')
        parts.append(text(x + 14, y + 138, "新增字段" if zh else "adds", 11, 600, SUB))
        for j, (slot, en, zh_text) in enumerate(st["adds"]):
            parts.append(text(x + 14, y + 158 + j * 26, "+ " + slot, 11.5, 700, BLUE[0]))
            parts.append(text(x + 26, y + 171 + j * 26, zh_text if zh else en, 10.5, 400, SUB))
        if i < len(STAGES) - 1:
            parts.append(poly_arrow([(x + card_w + 3, y + 27), (x + card_w + gap - 3, y + 27)], YELLOW[0]))
    y += 262 + 18
    # RCTD 的 R 旁路
    x0 = M + 1 * (card_w + gap)
    x1 = M + 2 * (card_w + gap) + card_w
    parts.append(rect(x0, y, x1 - x0, 62, fill=RED[1], stroke=RED[0], rx=12, sw=1.2, dash="5 4"))
    detour = ("RCTD 旁路（P1、P3）：02 的 counts → export_for_rctd.py → rctd/*.mtx → run_rctd.Rmd（R）→ weights.csv → import_rctd.py → 03" if zh else
              "RCTD detour (P1, P3): 02 counts → export_for_rctd.py → rctd/*.mtx → run_rctd.Rmd (R) → weights.csv → import_rctd.py → 03")
    for j, line in enumerate(wrap_text(detour, 12, x1 - x0 - 32)[:2]):
        parts.append(text(x0 + 16, y + 26 + j * 18, line, 12, 600 if j == 0 else 400, RED[0]))
    y += 62 + 34
    # 下游各读 05 的哪一部分
    parts.append(text(M + 6, y, "下游脚本读 05_niche.h5ad 的哪一部分" if zh else "Which part of 05_niche.h5ad each downstream script reads", 18, 700))
    y += 18
    box_h = len(CONSUMERS) * 74 - 12
    parts.append(rect(M, y, 230, box_h, fill=YELLOW[1], stroke=YELLOW[0], rx=16, sw=1.5))
    parts.append(icon("file", M + 115, y + box_h / 2 - 22, "acc", 24))
    parts.append(text(M + 115, y + box_h / 2 + 26, "05_niche.h5ad", 16, 700, YELLOW_DARK, "middle"))
    for i, (slot, en, zh_text, scripts) in enumerate(CONSUMERS):
        ry = y + i * 74
        parts.append(poly_arrow([(M + 230, y + box_h / 2), (M + 262, y + box_h / 2), (M + 262, ry + 31), (M + 286, ry + 31)], BLUE[0], 1.6))
        parts.append(rect(M + 290, ry, W - 2 * M - 290, 62, fill=BLUE[1], stroke=BLUE[0], rx=12, sw=1.2))
        parts.append(text(M + 306, ry + 26, slot, 13, 700, BLUE[0]))
        parts.append(text(M + 306, ry + 46, zh_text if zh else en, 11, 400, SUB))
        px = M + 640
        for name in scripts:
            label = ("全部图" if zh else "all figures") if name == "all figures" else name
            svg, pw = pill(px, ry + 19, label, "Py", 11)
            parts.append(svg); px += pw + 8
    y += box_h + 24
    note = ("layers['counts'] 只在质控和注释（含 R 的 RCTD）里用；下游分析都用 X（log 表达）。坐标一律用 obsm['spatial_um']（µm）。" if zh else
            "layers['counts'] is used only by QC and annotation (including RCTD in R); downstream analyses use X (log expression). Coordinates always come from obsm['spatial_um'] (µm).")
    parts.append(rect(M, y, W - 2 * M, 44, fill=GREY[1], stroke=BORDER, rx=12))
    parts.append(text(M + 18, y + 27, note, 12.5, 400, TEXT))
    return wrap(W, y + 64, "".join(parts))

# %%
# 步骤：写出全部 SVG
# 输入：PIPELINES
# 输出：docs/decision_{en,zh}.svg、docs/steps_{en,zh}.svg、docs/P*_{en,zh}.svg
for lang in ("en", "zh"):
    (OUT / f"decision_{lang}.svg").write_text(draw_decision(lang), encoding="utf-8")
    (OUT / f"steps_{lang}.svg").write_text(draw_steps(lang), encoding="utf-8")
    (OUT / f"dataflow_{lang}.svg").write_text(draw_dataflow(lang), encoding="utf-8")
    for p in PIPELINES:
        (OUT / f"{p['id']}_{lang}.svg").write_text(draw_pipeline(p, lang), encoding="utf-8")
print(sorted(f.name for f in OUT.glob("*.svg")))
