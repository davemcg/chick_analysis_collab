import scvelo as scv
import scanpy as sc
import argparse

parser = argparse.ArgumentParser(description='Run scVelo and save new adata object')
parser.add_argument('adata', help='Anndata object to run scVelo on')
parser.add_argument('adata_out', help='What to save the new adata object as')

parser.add_argument('--redo_pca', help='Redo HVG/PCA in anndata object', action='store_true')

args = parser.parse_args()
print(args)

adata = sc.read_h5ad(args.adata)

print(adata)

if args.redo_pca:
  print("\nRunning HVG and PCA\n")
  scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
  sc.tl.pca(adata)

sc.pp.neighbors(adata, n_pcs=30, n_neighbors=30)
# sc.tl.umap(adata)
scv.pp.moments(adata, n_pcs=None, n_neighbors=None)
print("Recover dynamics")
scv.tl.recover_dynamics(adata, n_jobs=10)

scv.tl.velocity(adata, mode="dynamical")
scv.tl.velocity_graph(adata, n_jobs=4)

adata.write_h5ad(args.adata_out)
