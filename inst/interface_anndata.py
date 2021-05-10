import anndata as ad
import scanpy as sp

def get_anndata(fname):
  adata = ad.read(fname)
  return(adata)
  
def cal_obsm(adata,obsm):
  if obsm == 'PCA':
    sp.tl.pca(adata,n_comps = 3)
  elif obsm == 'UMAP':
    try:
      adata.obsp['connectivities']
    except:
      sp.pp.neighbors(adata)
    sp.tl.umap(adata,n_components = 3)
  elif obsm == 'TSNE':
    sp.tl.tsne(adata)
  elif obsm == 'DRAW_GRAPH_FA':
    try:
      adata.obsp['connectivities']
    except:
      sp.pp.neighbors(adata)
    sp.tl.draw_graph(adata)
  elif obsm == 'DIFFMAP':
    try:
      adata.obsp['connectivities']
    except:
      sp.pp.neighbors(adata)
    sp.tl.diffmap(adata, n_comps = 3)
  return(adata)
  
def save_all_h5ad(adata, filename):
  adata.write(filename)


def save_sub_h5ad(adata, filename, key):
  sub_adata = adata[key]
  del(sub_adata.obsm)
  del(sub_adata.varm)
  del(sub_adata.uns)
  del(sub_adata.obsp)
  sub_adata.write(filename)

