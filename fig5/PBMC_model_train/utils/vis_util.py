import os
import numpy as np
import torch
import umap

import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from torch.autograd import Variable
import matplotlib.cm as cm
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import KNeighborsClassifier

def plot_pca_batch(anndat):
    batch = anndat.obs['batch'].unique().to_list()

    pca = PCA() #n_components=2)
    pca.fit(anndat.obsm['embed'])
    X_reduced = pca.transform(anndat.obsm['embed'])

    CAND_COLORS = cm.rainbow(np.linspace(0, 1, len(batch)))
    legend_elements = []
    for i in range(len(batch)):
        legend_elements.append(Line2D([i],[i],marker='o', color='w', label='%s'%(batch[i]), markerfacecolor=CAND_COLORS[i], markersize=6))

    plt.clf()
    fig, ax = plt.subplots(figsize=(5,5))

    scatter = ax.scatter(X_reduced[:, 0], X_reduced[:, 1], s=5, alpha=0.45, c=CAND_COLORS[anndat.obs['batch'].cat.codes])
    
    ax.set_xlabel("PCA1-%.2f%%"%(pca.explained_variance_ratio_[0]*100))
    ax.set_ylabel("PCA2-%.2f%%"%(pca.explained_variance_ratio_[1]*100))
    #ax.legend(handles=legend_elements, loc='upper left')
    ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1, 1))
    plt.show()



def plot_pca_label(anndat):
    labels = anndat.obs['cell_type'].unique().to_list()

    pca = PCA() #n_components=2)
    pca.fit(anndat.obsm['embed'])
    X_reduced = pca.transform(anndat.obsm['embed'])

    CAND_COLORS = cm.rainbow(np.linspace(0, 1, len(labels)))
    legend_elements = []
    for i in range(len(labels)):
        legend_elements.append(Line2D([i],[i],marker='o', color='w', label='%s'%(labels[i]), markerfacecolor=CAND_COLORS[i], markersize=6))

    plt.clf()
    fig, ax = plt.subplots(figsize=(5,5))
    
    scatter = ax.scatter(X_reduced[:, 0], X_reduced[:, 1], s=5, alpha=0.45, c=CAND_COLORS[anndat.obs['cell_type'].cat.codes])
    
    ax.set_xlabel("PCA1-%.2f%%"%(pca.explained_variance_ratio_[0]*100))
    ax.set_ylabel("PCA2-%.2f%%"%(pca.explained_variance_ratio_[1]*100))

    #ax.legend(handles=legend_elements, loc='upper left')
    ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1, 1))
    plt.show()


def umap_plot(representations, labels, min_dist=0.5, n_neigh=10):
    #pca = PCA(n_components=30)
    #pca.fit(representations[indeces])
    #X_reduced = pca.transform(representations[indeces])
    #X_reduced = np.multiply(X_reduced, pca.explained_variance_ratio_)

    embedding = umap.UMAP(n_neighbors=n_neigh,
                      min_dist=min_dist,
                      metric='euclidean')
    u = embedding.fit_transform(representations)

    CAND_COLORS = cm.rainbow(np.linspace(0, 1, 2))
    legend_elements = [Line2D([0], [0], marker='o', color='w', label='GEX',
                          markerfacecolor=CAND_COLORS[0], markersize=8),
                        Line2D([0], [0], marker='o', color='w', label='ATAC',
                          markerfacecolor=CAND_COLORS[1], markersize=8)]
    domains = np.zeros(len(u))
    domains[len(u)//2:] = 1
    domains = domains.astype('int')
    plt.clf()
    fig, ax = plt.subplots(figsize=(7,7))
    scatter = ax.scatter(u[:, 0], u[:, 1], s=5, alpha=0.6, c=CAND_COLORS[domains])
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    #legend1 = ax.legend(*scatter.legend_elements(),
    #                loc="upper left", title="Classes")
    #ax.add_artist(legend1)
    ax.legend(handles=legend_elements, loc='upper left')
    plt.show() #savefig('eval/%s/umap_domains_%s_%s.png'%(expname, lamb, margin))
    
    cell_type = np.unique(labels).tolist()
    nlabels = len(cell_type)

    CAND_COLORS = cm.rainbow(np.linspace(0, 1, nlabels))
    legend_elements = []
    for i in range(nlabels):
        legend_elements.append(Line2D([i],[i],marker='o', color='w', label='%s'%(cell_type[i]), markerfacecolor=CAND_COLORS[i], markersize=6))
        
    plt.clf()
    fig, ax = plt.subplots(figsize=(7,7))
    ax.scatter(u[:, 0], u[:, 1], s=5, alpha=0.6, c=CAND_COLORS[labels])
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1, 1))
    plt.show() #savefig('eval/%s/umap_labels_%s_%s.png'%(expname, lamb, margin))
    
    return u[:,0:2]

