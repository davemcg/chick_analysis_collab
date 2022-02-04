renv::init()
renv::install("reticulate")
renv::use_python() # 6 (~/git/chick_analysis_collab/renv/python/virtualenvs/renv-python-3.9/bin/python3)

pkgs <- c(
  "renv",
  "reticulate",
  "png",
  "ggplot2",
  "BiocManager",
  "Seurat",
  "tidyverse",
  "pals"
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
  "cellrank",
  'pandas==1.3.5'
)

reticulate::py_install(py_pkgs)

# only on github
remotes::install_github("mojaveazure/seurat-disk")
# need latest update that fixes import bug
remotes::install_github('theislab/zellkonverter')

renv::snapshot()
