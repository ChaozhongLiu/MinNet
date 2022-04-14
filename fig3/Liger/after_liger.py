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
meta = pd.read_csv('../bmmc.meta.s3s4.csv')
out_x = pd.read_csv('Liger_embed.rna.s3s4.%s.txt'%(randn),header=None).iloc[:,1:].to_numpy()
out_y = pd.read_csv('Liger_embed.atac.s3s4.%s.txt'%(randn),header=None).iloc[:,1:].to_numpy()
merge = np.concatenate((out_x, out_y), axis=0)
merge = pd.DataFrame(merge)
merge.insert(0, 'index', meta.index.to_numpy())
merge.to_csv('../results/raw/Liger_embed.s3s4.%s.txt'%(randn), index=False)


#-----------------------------------------------------------------------------------------
meta = pd.read_csv('../bmmc.meta.s4s3.csv')
out_x = pd.read_csv('Liger_embed.rna.s4s3.%s.txt'%(randn),header=None).iloc[:,1:].to_numpy()
out_y = pd.read_csv('Liger_embed.atac.s4s3.%s.txt'%(randn),header=None).iloc[:,1:].to_numpy()
merge = np.concatenate((out_x, out_y), axis=0)
merge = pd.DataFrame(merge)
merge.insert(0, 'index', meta.index.to_numpy())
merge.to_csv('../results/raw/Liger_embed.s4s3.%s.txt'%(randn), index=False)


#-----------------------------------------------------------------------------------------
meta = pd.read_csv('../bmmc.meta.s3s4.csv')
out_x = pd.read_csv('Liger_embed.rna.s3s4.%s.txt'%(randn),header=None).iloc[:,1:].to_numpy()
out_y = pd.read_csv('Liger_embed.atac.s3s4.%s.txt'%(randn),header=None).iloc[:,1:].to_numpy()

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

labels_x = meta['cell_type'].to_numpy()[0:out_x.shape[0]]
labels_y = meta['cell_type'].to_numpy()[out_x.shape[0]:]

pred_y = KNN_acc(out_x, out_y, labels_x, labels_y, nn=25)
pred_df = pd.DataFrame({'cell_type':labels_y,'prediction':pred_y})
pred_df.to_csv('../results/raw/Liger_prediction.s3s4.%s.txt'%(randn))

ct_df = pd.DataFrame({'cell_type':labels_y,
                     'prediction':pred_y})
n_ct = len(np.unique(labels_y))
ct_heatmap = pd.DataFrame(np.zeros((n_ct,n_ct)))
ct_heatmap.index = meta['cell_type'][out_x.shape[0]:].unique().tolist()
ct_heatmap.columns = meta['cell_type'][out_x.shape[0]:].unique().tolist()
for i in range(n_ct):
    for j in range(n_ct):
        tmp = ct_df.loc[(ct_df['cell_type'] == ct_heatmap.index[i]),:]
        ct_heatmap.iloc[i,j] = np.mean(tmp['prediction'] == ct_heatmap.columns[j])

ct_heatmap.index = ct_heatmap.index + ' - ' + ct_df['cell_type'].value_counts().astype(str)[ct_heatmap.index]

acc_df = pd.DataFrame({'cell_type':ct_heatmap.columns.tolist(),
                       'GLUE': np.diag(ct_heatmap)}) #wrong name
acc_df.to_csv('../results/raw/Liger_transfer_acc.s3s4.%s.csv'%(randn))



#-----------------------------------------------------------------------------------------
meta = pd.read_csv('../bmmc.meta.s4s3.csv')
out_x = pd.read_csv('Liger_embed.rna.s4s3.%s.txt'%(randn),header=None).iloc[:,1:].to_numpy()
out_y = pd.read_csv('Liger_embed.atac.s4s3.%s.txt'%(randn),header=None).iloc[:,1:].to_numpy()

labels_x = meta['cell_type'].to_numpy()[0:out_x.shape[0]]
labels_y = meta['cell_type'].to_numpy()[out_x.shape[0]:]

pred_y = KNN_acc(out_x, out_y, labels_x, labels_y, nn=25)
pred_df = pd.DataFrame({'cell_type':labels_y,'prediction':pred_y})
pred_df.to_csv('../results/raw/Liger_prediction.s4s3.%s.txt'%(randn))

ct_df = pd.DataFrame({'cell_type':labels_y,
                     'prediction':pred_y})
n_ct = len(np.unique(labels_y))
ct_heatmap = pd.DataFrame(np.zeros((n_ct,n_ct)))
ct_heatmap.index = meta['cell_type'][out_x.shape[0]:].unique().tolist()
ct_heatmap.columns = meta['cell_type'][out_x.shape[0]:].unique().tolist()
for i in range(n_ct):
    for j in range(n_ct):
        tmp = ct_df.loc[(ct_df['cell_type'] == ct_heatmap.index[i]),:]
        ct_heatmap.iloc[i,j] = np.mean(tmp['prediction'] == ct_heatmap.columns[j])

ct_heatmap.index = ct_heatmap.index + ' - ' + ct_df['cell_type'].value_counts().astype(str)[ct_heatmap.index]

acc_df = pd.DataFrame({'cell_type':ct_heatmap.columns.tolist(),
                       'GLUE': np.diag(ct_heatmap)}) #wrong name
acc_df.to_csv('../results/raw/Liger_transfer_acc.s4s3.%s.csv'%(randn))










