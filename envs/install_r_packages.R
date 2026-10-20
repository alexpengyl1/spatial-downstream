# R 环境：RCTD（spacexr）、BANKSY，以及读写交换文件需要的包
# 已验证：macOS arm64，R 4.6.1，Bioconductor 3.23，spacexr 2.2.1，Banksy 1.8.1
# 为什么不用 conda：bioconda 上的 r-spacexr 是给 R 4.3 编译的，bioconductor-banksy 是给 R 4.5 编译的，两者装不进同一个环境；
#                  直接装到系统 R 最省事（流程脚本里用 R_CMD="Rscript"）
# 运行：`Rscript install_r_packages.R`
# R 步骤都是 R Markdown 报告（.Rmd）：渲染需要 pandoc，RStudio / Quarto 自带；命令行找不到时设置 RSTUDIO_PANDOC 指向它

# 步骤：CRAN 包（包括 spacexr 的依赖：install_github 不一定能自动装全）
# 输入：无
# 输出：yaml、Matrix、remotes、BiocManager 和 spacexr 的依赖
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 4)
install.packages(c("yaml", "Matrix", "remotes", "BiocManager", "rmarkdown", "knitr",
                   "quadprog", "doParallel", "foreach", "fields",
                   "pals", "reshape2", "CompQuadForm", "Rfast", "locfdr", "metafor"))

# 步骤：BANKSY 及依赖（Bioconductor 有 macOS 预编译包）
# 输入：无
# 输出：Banksy、SpatialExperiment、scuttle、scran
BiocManager::install(c("Banksy", "SpatialExperiment", "scuttle", "scran"), ask = FALSE, update = FALSE)

# 步骤：RCTD
# 输入：无
# 输出：spacexr
# 如果报 "download from 'https://api.github.com/...' failed"（git 凭据里存的 GitHub PAT 失效时常见），改在终端运行：
#   git clone --depth 1 https://github.com/dmcable/spacexr.git && R CMD INSTALL --no-build-vignettes spacexr
remotes::install_github("dmcable/spacexr", build_vignettes = FALSE, upgrade = "never")

# 步骤：确认安装成功
# 输入：无
# 输出：打印的版本号
packageVersion("spacexr")
packageVersion("Banksy")

# 步骤：可选，在 R 里直接读写 h5ad（见 references/r_python_bridge.md 的方式二）
# 输入：无
# 输出：zellkonverter
# BiocManager::install("zellkonverter")

# 步骤：可选，seurat_to_mtx.Rmd 需要（读 Seurat 对象）
# 输入：无
# 输出：SeuratObject
# install.packages("SeuratObject")
