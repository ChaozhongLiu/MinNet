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
#neighborhood
meta = pd.read_csv('../cite.rna.meta.csv')
out_x = pd.read_csv('Liger_embed.rna.%s.txt'%(randn),header=None)
out_y = pd.read_csv('Liger_embed.adt.%s.txt'%(randn),header=None)

match_matrix = pairwise_distances(out_x.iloc[:,1:].to_numpy(), 
                                  out_y.iloc[:,1:].to_numpy(), metric="euclidean")
#match_matrix = np.exp(-match_matrix)

batch_code = meta['batch'].to_numpy()

match_matrix_1 = match_matrix[batch_code=='s4d1',:]
match_matrix_1 = match_matrix_1[:,batch_code=='s4d1']

match_matrix_2 = match_matrix[batch_code=='s4d8',:]
match_matrix_2 = match_matrix_2[:,batch_code=='s4d8']

match_matrix_3 = match_matrix[batch_code=='s4d9',:]
match_matrix_3 = match_matrix_3[:,batch_code=='s4d9']

neighbor_1 = np.argsort(match_matrix_1, axis=1, kind='stable')
neighbor_1 = np.argsort(neighbor_1, axis=1, kind='stable')

neighbor_1 = np.argsort(match_matrix_1, axis=1, kind='stable')
neighbor_1 = np.argsort(neighbor_1, axis=1, kind='stable')
writer = open('../results/raw/neighborhood_s4d1_Liger.%s.txt'%(randn),'w')
writer.write('\n'.join(np.diag(neighbor_1).astype('str').tolist()))
writer.close()


neighbor_2 = np.argsort(match_matrix_2, axis=1, kind='stable')
neighbor_2 = np.argsort(neighbor_2, axis=1, kind='stable')
writer = open('../results/raw/neighborhood_s4d8_Liger.%s.txt'%(randn),'w')
writer.write('\n'.join(np.diag(neighbor_2).astype('str').tolist()))
writer.close()


neighbor_3 = np.argsort(match_matrix_3, axis=1, kind='stable')
neighbor_3 = np.argsort(neighbor_3, axis=1, kind='stable')
writer = open('../results/raw/neighborhood_s4d9_Liger.%s.txt'%(randn),'w')
writer.write('\n'.join(np.diag(neighbor_3).astype('str').tolist()))
writer.close()



#-----------------------------------------------------------------------------------------
#label transfer
meta = pd.read_csv('../cite.rna.meta.csv')
out_x = pd.read_csv('Liger_embed.adt.%s.txt'%(randn),header=None).iloc[:,1:].to_numpy()
out_y = pd.read_csv('Liger_embed.rna.%s.txt'%(randn),header=None).iloc[:,1:].to_numpy()

input_mod1 = ad.read_h5ad('../data/Cite_GEX.test.h5ad')


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


merge = np.concatenate((out_x, out_y), axis=0)
labels = meta['cell_type'].astype('category').cat.codes.to_numpy()

pred_y = KNN_acc(merge, labels, nn=200)
input_mod1.obs['cell_type_name'] = input_mod1.obs['cell_type']
input_mod1.obs['cell_type'] = input_mod1.obs['cell_type'].cat.codes
cell_type_code = input_mod1.obs[['cell_type','cell_type_name']]
cell_type_code = pd.DataFrame({'code':cell_type_code['cell_type'].unique(),
                              'name':cell_type_code['cell_type_name'].unique()})
cell_type_code = cell_type_code.sort_values('code')
cell_type_code.index = cell_type_code['code']
pred_y = cell_type_code.loc[pred_y,'name'].to_numpy()
labels_x = input_mod1.obs['cell_type_name'].to_numpy()

pred_df = pd.DataFrame({'cell_type':labels_x,'prediction':pred_y})
pred_df.to_csv('../results/raw/Liger_prediction.%s.txt'%(randn))

ct_df = pd.DataFrame({'cell_type':labels_x,
                     'prediction':pred_y})

ct_heatmap = pd.DataFrame(np.zeros((23,23)))
ct_heatmap.index = input_mod1.obs['cell_type_name'].unique().tolist()
ct_heatmap.columns = input_mod1.obs['cell_type_name'].unique().tolist()
for i in range(23):
    for j in range(23):
        tmp = ct_df.loc[(ct_df['cell_type'] == ct_heatmap.index[i]),:]
        ct_heatmap.iloc[i,j] = np.mean(tmp['prediction'] == ct_heatmap.columns[j])
        
ct_heatmap.index = ct_heatmap.index + ' - ' + ct_df['cell_type'].value_counts().astype(str)[ct_heatmap.index]
mask = ct_heatmap.isnull()

acc_df = pd.DataFrame({'cell_type':ct_heatmap.columns.tolist(),
                       'Liger': np.diag(ct_heatmap)})
acc_df.to_csv('../results/raw/Liger_transfer_acc.%s.csv'%(randn))




#-----------------------------------------------------------------------------------------
meta = pd.read_csv('../cite.rna.meta.csv')
out_x = pd.read_csv('Liger_embed.rna.%s.txt'%(randn),header=None).iloc[:,1:].to_numpy()
out_y = pd.read_csv('Liger_embed.adt.%s.txt'%(randn),header=None).iloc[:,1:].to_numpy()
merge = np.concatenate((out_x, out_y), axis=0)
merge = pd.DataFrame(merge)
merge.insert(0, 'index', np.concatenate((meta.index.to_numpy(), 
                                         meta.index.to_numpy()), axis=0))
merge.to_csv('../results/raw/Liger_embed.%s.txt'%(randn), index=False)




