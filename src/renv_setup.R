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
  "SingleCellExperiment",
  "scater",
  "multtest"
)

renv::install(pkgs)


py_pkgs <- c(
  "scanpy",
  "python-igraph",
  "louvain",
  "scvelo",
  "cellrank"
)

reticulate::py_install(py_pkgs)

renv::snapshot()
