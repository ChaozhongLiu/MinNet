import anndata
import networkx as nx
import scanpy as sc
import scglue
from matplotlib import rcParams
import numpy as np
import scipy as sp
import seaborn as sns
import itertools
import pandas as pd
import sys

randn = int(sys.argv[1])

rna = anndata.read_h5ad("rna_preprocessed.h5ad")
atac = anndata.read_h5ad("atac_preprocessed.h5ad")
graph = nx.read_graphml("prior.graphml.gz")


scglue.models.configure_dataset(
    rna, "NB", use_highly_variable=True,
    use_layer="raw", use_rep="X_pca"
)

scglue.models.configure_dataset(
    atac, "NB", use_highly_variable=True,
    use_rep="X_lsi"
)

graph = graph.subgraph(itertools.chain(
    rna.var.query("highly_variable").index,
    atac.var.query("highly_variable").index
))

glue = scglue.models.SCGLUEModel(
    {"rna": rna, "atac": atac}, sorted(graph.nodes),
    random_seed=randn
)

glue.compile()

glue.fit(
    {"rna": rna, "atac": atac},
    graph, edge_weight="weight", edge_sign="sign",
    directory="glue_s4s3"
)

glue.save("glue_s4s3/final.%s.dill"%(randn))



rna.obsm["X_glue"] = glue.encode_data("rna", rna)
atac.obsm["X_glue"] = glue.encode_data("atac", atac)
combined = anndata.AnnData(
    obs=pd.concat([rna.obs, atac.obs], join="inner"),
    obsm={"X_glue": np.concatenate([rna.obsm["X_glue"], atac.obsm["X_glue"]])}
)

combined.obs['domain'] = 'RNA'
combined.obs['domain'][rna.n_obs:] = 'ATAC'

sc.pp.neighbors(combined, use_rep="X_glue", metric="cosine")
sc.tl.umap(combined)
#sc.pl.umap(combined, color=["cell_type", "domain", "batch"], wspace=0.65)

merge = pd.DataFrame(combined.obsm['X_glue'])
merge.insert(0, 'index', combined.obs.index.to_numpy())
merge.to_csv('../results/raw/GLUE_embed.s4s3.%s.txt'%(randn), index=False)



from sklearn.decomposition import PCA
from sklearn.neighbors import KNeighborsClassifier

def KNN_acc(out_x, out_y, labels_x, labels_y, nn=5):
    neigh = KNeighborsClassifier(n_neighbors=nn, weights='distance')
    neigh.fit(out_x, labels_x)
    labels_pred = neigh.predict(out_y)
    #print(check)
    acc = np.mean(labels_y == labels_pred)
    print("Accuracy is %s"%(acc))
    return labels_pred

out_x = merge.iloc[0:rna.shape[0],1:].to_numpy()
out_y = merge.iloc[rna.shape[0]:,1:].to_numpy()
labels_x = rna.obs['cell_type'].to_numpy()
labels_y = atac.obs['cell_type'].to_numpy()

pred_y = KNN_acc(out_x, out_y, labels_x, labels_y, nn=50)

pred_df = pd.DataFrame({'cell_type':labels_y,'prediction':pred_y})
pred_df.to_csv('../results/raw/GLUE_prediction.s4s3.%s.txt'%(randn))

ct_df = pd.DataFrame({'cell_type':labels_y,
                     'prediction':pred_y})
n_ct = len(np.unique(labels_y))
ct_heatmap = pd.DataFrame(np.zeros((n_ct,n_ct)))
ct_heatmap.index = atac.obs['cell_type'].unique().tolist()
ct_heatmap.columns = atac.obs['cell_type'].unique().tolist()
for i in range(n_ct):
    for j in range(n_ct):
        tmp = ct_df.loc[(ct_df['cell_type'] == ct_heatmap.index[i]),:]
        ct_heatmap.iloc[i,j] = np.mean(tmp['prediction'] == ct_heatmap.columns[j])

ct_heatmap.index = ct_heatmap.index + ' - ' + ct_df['cell_type'].value_counts().astype(str)[ct_heatmap.index]

acc_df = pd.DataFrame({'cell_type':ct_heatmap.columns.tolist(),
                       'GLUE': np.diag(ct_heatmap)}) #wrong name
acc_df.to_csv('../results/raw/GLUE_transfer_acc.s4s3.%s.csv'%(randn))
        
        

