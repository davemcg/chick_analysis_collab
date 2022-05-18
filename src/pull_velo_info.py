import scvelo as scv
import scanpy as sc
import cellrank as cr

def pull_velo_info(adata_path, output_prefix):
  # adata_path = '/Users/mcgaugheyd/data/chick_miruna/m006_scVelo.h5ad'
  adata = sc.read_h5ad(adata_path)
  velo = scv.DataFrame(adata.layers['velocity'])
  velo.to_csv(output_prefix + '.velo.csv')
  genes = adata.var['velocity_genes']
  genes.to_csv(output_prefix + '.velogenes.csv')
  obs = scv.DataFrame(adata.obs)
  obs.to_csv(output_prefix + '.obs.csv')

pull_velo_info('/Users/mcgaugheyd/data/chick_miruna/m006_scVelo.h5ad', '/Users/mcgaugheyd/git/chick_analysis_collab/data/m006_scVelo')
pull_velo_info('/Users/mcgaugheyd/data/chick_miruna/m007_scVelo.h5ad', '/Users/mcgaugheyd/git/chick_analysis_collab/data/m007_scVelo')
pull_velo_info('/Users/mcgaugheyd/data/chick_miruna/m008_scVelo.h5ad', '/Users/mcgaugheyd/git/chick_analysis_collab/data/m008_scVelo')
  
