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
    directory="glue"
)

glue.save("glue/final.%s.dill"%(randn))

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
merge.to_csv('../results/raw/GLUE_embed.%s.txt'%(randn), index=False)


#------------------------------------------------------------------------------------------
from sklearn.preprocessing import normalize
from sklearn.metrics import pairwise_distances

match_matrix = pairwise_distances(rna.obsm['X_glue'], 
                                  atac.obsm['X_glue'])
#match_matrix = np.exp(-match_matrix)

batch_code = rna.obs['batch'].to_numpy()

match_matrix_1 = match_matrix[batch_code=='s1d2',:]
match_matrix_1 = match_matrix_1[:,batch_code=='s1d2']

match_matrix_2 = match_matrix[batch_code=='s3d7',:]
match_matrix_2 = match_matrix_2[:,batch_code=='s3d7']

neighbor_1 = np.argsort(match_matrix_1, axis=1, kind='stable')
neighbor_1 = np.argsort(neighbor_1, axis=1, kind='stable')

writer = open('../results/raw/neighborhood_s1d2_glue.%s.txt'%(randn),'w')
writer.write('\n'.join(np.diag(neighbor_1).astype('str').tolist()))
writer.close()

neighbor_2 = np.argsort(match_matrix_2, axis=1, kind='stable')
neighbor_2 = np.argsort(neighbor_2, axis=1, kind='stable')

writer = open('../results/raw/neighborhood_s3d7_glue.%s.txt'%(randn),'w')
writer.write('\n'.join(np.diag(neighbor_2).astype('str').tolist()))
writer.close()


#------------------------------------------------------------------------------------------
from sklearn.decomposition import PCA
from sklearn.neighbors import KNeighborsClassifier

def KNN_acc(representations, labels, nn=5):
    indeces = np.arange(representations.shape[0])
    #representations = StandardScaler().fit_transform(representations)
    #pca = PCA(n_components=10)
    #pca.fit(representations[indeces])
    #X_reduced = pca.transform(representations[indeces])
    X_reduced = representations
    neigh = KNeighborsClassifier(n_neighbors=nn, weights='distance')
    n = int(len(representations) // 2)
    neigh.fit(X_reduced[0:n,:], labels)
    ylabel = neigh.predict(X_reduced[n:,:])
    check = len(ylabel) == n
    #print(check)
    acc = np.mean(labels == ylabel)
    print("Accuracy is %s"%(acc))
    return ylabel


labels = rna.obs['cell_type'].cat.codes.to_numpy()
merge = np.concatenate((rna.obsm['X_glue'], atac.obsm['X_glue']), axis=0)

pred_y = KNN_acc(merge, labels, nn=25)

rna.obs['cell_type_name'] = rna.obs['cell_type']
rna.obs['cell_type'] = rna.obs['cell_type'].cat.codes
cell_type_code = rna.obs[['cell_type','cell_type_name']]
cell_type_code = pd.DataFrame({'code':cell_type_code['cell_type'].unique(),
                              'name':cell_type_code['cell_type_name'].unique()})
cell_type_code = cell_type_code.sort_values('code')
cell_type_code.index = cell_type_code['code']

pred_y = cell_type_code.loc[pred_y,'name'].to_numpy()
labels_x = rna.obs['cell_type_name'].to_numpy()

pred_df = pd.DataFrame({'cell_type':labels_x,'prediction':pred_y})
pred_df.to_csv('../results/raw/GLUE_prediction.%s.txt'%(randn))

ct_df = pd.DataFrame({'cell_type':labels_x,
                     'prediction':pred_y})

ct_heatmap = pd.DataFrame(np.zeros((22,22)))
ct_heatmap.index = rna.obs['cell_type_name'].unique().tolist()
ct_heatmap.columns = rna.obs['cell_type_name'].unique().tolist()
for i in range(22):
    for j in range(22):
        tmp = ct_df.loc[(ct_df['cell_type'] == ct_heatmap.index[i]),:]
        ct_heatmap.iloc[i,j] = np.mean(tmp['prediction'] == ct_heatmap.columns[j])
        
ct_heatmap.index = ct_heatmap.index + ' - ' + ct_df['cell_type'].value_counts().astype(str)[ct_heatmap.index]
sns.set(rc={"figure.dpi":200, 'savefig.dpi':100},font_scale=0.4)
mask = ct_heatmap.isnull()
#sns.heatmap(ct_heatmap, mask=mask, cmap="BuPu")

acc_df = pd.DataFrame({'cell_type':ct_heatmap.columns.tolist(),
                       'glue': np.diag(ct_heatmap)}) #wrong name

acc_df.to_csv('../results/raw/GLUE_transfer_acc.%s.csv'%(randn))





