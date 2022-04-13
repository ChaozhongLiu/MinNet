import anndata as ad
import numpy as np
import scipy as sp
import pandas as pd
import seaborn as sns
import sys
import torch
import os
import matplotlib.pyplot as plt

from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import normalize
from sklearn.metrics import pairwise_distances

randn = int(sys.argv[1])


#-----------------------------------------------------------------------------------------
meta = pd.read_csv('../data/rna.meta.csv',sep='\t')
out_x = pd.read_csv('Liger_embed.rna.%s.txt'%(randn),header=None)
out_y = pd.read_csv('Liger_embed.atac.%s.txt'%(randn),header=None)

match_matrix = pairwise_distances(out_x.iloc[:,1:].to_numpy(), 
                                  out_y.iloc[:,1:].to_numpy(), metric="euclidean")

neighbor_1 = np.argsort(match_matrix, axis=1, kind='stable')
neighbor_1 = np.argsort(neighbor_1, axis=1, kind='stable')

writer = open('../results/raw/neighborhood_Liger.%s.txt'%(randn),'w')
writer.write('\n'.join(np.diag(neighbor_1).astype('str').tolist()))
writer.close()



#-----------------------------------------------------------------------------------------
meta = pd.read_csv('../data/rna.meta.csv',sep='\t')
out_x = pd.read_csv('Liger_embed.rna.%s.txt'%(randn),header=None).iloc[:,1:].to_numpy()
out_y = pd.read_csv('Liger_embed.atac.%s.txt'%(randn),header=None).iloc[:,1:].to_numpy()

input_mod1 = ad.read_h5ad('../data/PBMC_10X_GEX.h5ad')

from sklearn.decomposition import PCA
from sklearn.neighbors import KNeighborsClassifier

def KNN_acc(representations, labels, nn=5):
    indeces = np.arange(representations.shape[0])
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

merge = np.concatenate((out_x, out_y), axis=0)
labels = meta['seurat_annotations'].to_numpy()

pred_y = KNN_acc(merge, labels, nn=40)

input_mod1.obs['cell_type_name'] = input_mod1.obs['seurat_annotations']
input_mod1.obs['cell_type'] = input_mod1.obs['seurat_annotations'].cat.codes
cell_type_code = input_mod1.obs[['cell_type','cell_type_name']]
cell_type_code = pd.DataFrame({'code':cell_type_code['cell_type'].unique(),
                              'name':cell_type_code['cell_type_name'].unique()})
cell_type_code = cell_type_code.sort_values('code')
cell_type_code.index = cell_type_code['code']

pred_df = pd.DataFrame({'cell_type':labels,'prediction':pred_y})
pred_df.to_csv('../results/raw/Liger_prediction.%s.txt'%(randn))

ct_df = pd.DataFrame({'cell_type':labels,
                     'prediction':pred_y})

ct_heatmap = pd.DataFrame(np.zeros((19,19)))
ct_heatmap.index = input_mod1.obs['cell_type_name'].unique().tolist()
ct_heatmap.columns = input_mod1.obs['cell_type_name'].unique().tolist()
for i in range(19):
    for j in range(19):
        tmp = ct_df.loc[(ct_df['cell_type'] == ct_heatmap.index[i]),:]
        ct_heatmap.iloc[i,j] = np.mean(tmp['prediction'] == ct_heatmap.columns[j])
        
        
ct_heatmap.index = ct_heatmap.index + ' - ' + ct_df['cell_type'].value_counts().astype(str)[ct_heatmap.index]
sns.set(rc={"figure.dpi":200, 'savefig.dpi':100},font_scale=0.4)
mask = ct_heatmap.isnull()
sns.heatmap(ct_heatmap, mask=mask, cmap="BuPu")

acc_df = pd.DataFrame({'cell_type':ct_heatmap.columns.tolist(),
                       'Liger': np.diag(ct_heatmap)})

acc_df.to_csv('../results/raw/Liger_transfer_acc.%s.csv'%(randn))

#-----------------------------------------------------------------------------------------
meta = pd.read_csv('../data/rna.meta.csv',sep='\t')
out_x = pd.read_csv('Liger_embed.rna.%s.txt'%(randn),header=None).iloc[:,1:].to_numpy()
out_y = pd.read_csv('Liger_embed.atac.%s.txt'%(randn),header=None).iloc[:,1:].to_numpy()

merge = np.concatenate((out_x, out_y), axis=0)
merge = pd.DataFrame(merge)
merge.insert(0, 'index', np.concatenate((meta.index.to_numpy(), 
                                         meta.index.to_numpy()), axis=0))
merge.to_csv('../results/raw/Liger_embed.%s.txt'%(randn), index=False)



#-----------------------------------------------------------------------------------------
meta = pd.read_csv('../data/rna.meta.csv',sep='\t')
out_x = pd.read_csv('Liger_UMAPembed.rna.%s.txt'%(randn),header=None).iloc[:,1:].to_numpy()
out_y = pd.read_csv('Liger_UMAPembed.atac.%s.txt'%(randn),header=None).iloc[:,1:].to_numpy()
merge = np.concatenate((out_x, out_y), axis=0)
merge = pd.DataFrame(merge)
merge.insert(0, 'index', np.concatenate((meta.index.to_numpy(), 
                                         meta.index.to_numpy()), axis=0))
merge.to_csv('../results/raw/Liger_UMAP.%s.txt'%(randn), index=False)






