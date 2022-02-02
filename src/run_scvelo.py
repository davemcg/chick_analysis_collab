import scvelo as scv
import scanpy as sc
import argparse

parser = argparse.ArgumentParser(description='Run scVelo and save new adata object')
parser.add_argument('--use_existing_pca', help='Use exising PCA, (Y)es or (N)o')

adata = sc.read_h5ad('/Users/mcgaugheyd/data/chick_miruna/m007.h5ad')

if args[1] != 'use_exisiting_pca':
  scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
  sc.tl.pca(adata)

sc.pp.neighbors(adata, n_pcs=30, n_neighbors=30)
# sc.tl.umap(adata)
scv.pp.moments(adata, n_pcs=None, n_neighbors=None)

scv.tl.recover_dynamics(adata, n_jobs=10)

scv.tl.velocity(adata, mode="dynamical", n_jobs=4)
scv.tl.velocity_graph(adata)

adata.write_h5ad('/Users/mcgaugheyd/data/chick_miruna/m007_scvelo.h5ad')
