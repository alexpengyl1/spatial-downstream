# 原始示例脚本（第 1–28 课）→ 本 skill 脚本

| 课 | 课程脚本 | 本 skill | 主要改动 |
|---|---|---|---|
| 1 | scanpy_multi.py | B_preprocess/preprocess_rna.py | 去掉 argparse 和 regress_out；保留全部基因 |
| 2 | run_infercnv.R、Spatial.CNV.py | G_special/cnv_infercnvpy.py | 改用 infercnvpy；去掉 Jupyter 魔法命令 |
| 2 | cytotrace.r、mutinichenet.r、scvelo.py | 未收录（单细胞分析，见 r_python_bridge.md） | — |
| 3 | scanpy.spatial.py | A_load/load_visium.py、load_stereo.py + preprocess_rna.py | 修复未定义的 plt / outdir / path |
| 4 | RCTD_visium.r | C_annotate/export_for_rctd.py → run_rctd.Rmd → import_rctd.py | 原脚本保存的是 spatial 而不是 RCTD 结果；改为文件交换 |
| 4 | Spatial_cell2location.py | C_annotate/run_cell2location.py | 修复参数名拼写、缩进、未导入的 ad |
| 5 | RCTD.CC.r、Spatial.niche.py | D_domain_niche/niche_composition.py（spot 分支） | 统一成 KMeans；修复 np / outdir |
| 6 | spatial.colocalization.py | E_neighborhood/misty_colocalization.py | bandwidth 改为 µm |
| 7 | COMMOT.py | F_communication/commot_flow.py | 修复 args / argv；距离改为 µm；不再在循环里覆盖同一个图 |
| 7 | SpaGene.r | G_special/svg_moran.py + F_communication/lr_spatial_liana.py | Python 替代 |
| 8 | spatial.cnv.py | G_special/cnv_infercnvpy.py | 增加 GTF 基因位置；用 cnv_score |
| 9 | hopspot.py、spatial.hotspot.py | G_special/hotspot_getis.py | 改用 esda 的 Getis-Ord；距离用 KD 树 |
| 10 | spatial.tendency.py | G_special/spatial_tendency_soapy.py | 修复缩进；单样本、µm |
| 11 | scanpy.spatial.HD.py、tissue_positions.py | A_load/load_visium_hd.py | 直接读 parquet 坐标；去掉写死的样本名 |
| 12 | banksy.HD.py | D_domain_niche/domain_cellcharter.py（Py）或 domain_banksy.Rmd（R） | Python 版 BANKSY 依赖本地路径，改用 CellCharter 或 Bioconductor 版 |
| 13 | RCTD.CC.r | 同第4课，rctd_mode: doublet | 去掉写死的元数据路径 |
| 14 | scniche.hd.py | D_domain_niche/niche_composition.py（cell 分支） | 全距离矩阵改为稀疏邻居图；修复列名不一致 |
| 15 | coor.hd.py | E_neighborhood/nhood_enrichment.py | 加上 interaction_matrix 和 co_occurrence |
| 16 | coor.hd.diff.py | E_neighborhood/diff_nhood.py | 改为以样本为重复单位：每个样本做邻域富集，两组 z 分数做 Welch t 检验 |
| 17 | squidpy_cci.py | F_communication/lr_spatial_liana.py | ligrec 不使用空间信息，改用 bivariate |
| 18 | spatial.nmf.py | G_special/nmf_programs.py | 输出 top 基因表；先下采样拟合 |
| 19 | squidpy.foundation.py | A_load/load_xenium.py、load_cosmx.py + preprocess_rna.py | 修复未定义的 xenium()；CosMx 用 global 坐标 |
| 20 | label_transfer.py | C_annotate/label_transfer.py | 手写 KNN 投票改用 sklearn；修复传入字符串的问题 |
| 21 | scimap.cn.py | D_domain_niche/niche_composition.py | 不依赖 scimap，结果相同（邻域计数 + KMeans） |
| 22 | squidpy.Delaunay.py、cellcharter.neighbor.py | E_neighborhood/nhood_enrichment.py | 对角线置空代替逐行归一化 |
| 23 | cellcharter.diff.py | E_neighborhood/diff_nhood.py | CellCharter 0.3.5 的 diff_nhood_enrichment 在 n_jobs=1 时报错、多进程时卡住，改用同一脚本的样本级检验 |
| 24 | scimap.subtype.NMF.py | G_special/nmf_subtype.py | 增加亚型邻域组成图和差异基因表 |
| 25 | squidpy.neighbor.lr.py | F_communication/ligrec_neighborhood.py | 稀疏矩阵取邻居，不用 Python 循环；squidpy 1.6.5 的 ligrec p 值全为 NaN，改用 LIANA+ 的 CellPhoneDB |
| 26 | codex.py | A_load/load_codex.py + B_preprocess/preprocess_protein.py | 增加 arcsinh 转换 |
| 27 | SOAPy.pro.py、pro.RCN.py | D_domain_niche/niche_composition.py | 修复循环内重复赋值 |
| 28 | （只有 PDF） | references/method_matrix.md、units_and_params.md | 邻域定义方法的对比 |
