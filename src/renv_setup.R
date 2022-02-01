renv::init()
renv::install("reticulate")
renv::use_python() # 2 (/usr/local/Cellar/python@3.9/3.9.4/bin/python3.9)

pkgs <- c(
  "renv",
  "reticulate",
  "png",
  "ggplot2",
  "BiocManager",
  "Seurat"
)

bioc_pkgs <- c(
  "bioc::SingleCellExperiment",
  "bioc::scater",
  "bioc::multtest"
)

renv::install(pkgs)
renv::install(bioc_pkgs)


py_pkgs <- c(
  "scanpy",
  "python-igraph",
  "louvain",
  "scvelo",
  "cellrank"
)

reticulate::py_install(py_pkgs)

renv::snapshot()
remotes::install_github("mojaveazure/seurat-disk")
