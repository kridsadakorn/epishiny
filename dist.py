import anndata as ad
import scanpy as sp

def get_anndata(fname):
  adata = ad.read(fname)
  return(adata)
  
def cal_obsm(adata,obsm):
  if obsm == 'PCA':
    sp.tl.pca(adata,n_comps = 3)
  elif obsm == 'UMAP':
    sp.tl.umap(adata,n_components = 3)
  elif obsm == 'TSNE':
    sp.tl.tsne(adata)
  elif obsm == 'DRAW_GRAPH_FA':
    sp.tl.draw_graph(adata)
  return(adata)

