import os
import numpy as np
import torch
import umap
import umap.plot

import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from torch.autograd import Variable
import matplotlib.cm as cm
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import KNeighborsClassifier

def plot_pca(representations, labels, ct_dic, batch, batch_dic):
    pca = PCA() #n_components=2)
    pca.fit(representations)
    u = pca.transform(representations)
    
    #----------------------------------------------------------
    cell_type = np.unique(labels).tolist()
    nlabels = len(cell_type)
    CAND_COLORS = cm.rainbow(np.linspace(0, 1, nlabels))
    legend_elements = []
    for i in range(nlabels):
        legend_elements.append(Line2D([i],[i],marker='o', color='w', label='%s'%(ct_dic[i]), markerfacecolor=CAND_COLORS[i], markersize=6))
    
    plt.clf()
    fig, ax = plt.subplots(figsize=(5,5))
    scatter = ax.scatter(u[:, 0], u[:, 1], s=3, alpha=0.5, c=CAND_COLORS[labels])
    ax.set_xlabel("PCA1-%.2f%%"%(pca.explained_variance_ratio_[0]*100))
    ax.set_ylabel("PCA2-%.2f%%"%(pca.explained_variance_ratio_[1]*100))
    #legend1 = ax.legend(*scatter.legend_elements(),
    #                loc="upper left", title="Classes")
    #ax.add_artist(legend1)
    ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1, 1))
    plt.show()
    
    
    #------------------------------------------------------------------------------
    batch_list = np.unique(batch).tolist()
    nlabels = len(batch_list)
    CAND_COLORS = cm.rainbow(np.linspace(0, 1, nlabels+2))
    legend_elements = []
    for i in range(nlabels):
        legend_elements.append(Line2D([i],[i],marker='o', color='w', label='%s'%(batch_dic[i]), markerfacecolor=CAND_COLORS[i+1], markersize=6))
    
    plt.clf()
    fig, ax = plt.subplots(figsize=(5,5))
    ax.scatter(u[:, 0], u[:, 1], s=3, alpha=0.5, c=CAND_COLORS[batch+1])
    ax.set_xlabel("PCA1-%.2f%%"%(pca.explained_variance_ratio_[0]*100))
    ax.set_ylabel("PCA2-%.2f%%"%(pca.explained_variance_ratio_[1]*100))
    ax.legend(handles=legend_elements, loc='upper left')
    plt.show()


def umap_plot(representations, labels, ct_dic, batch, batch_dic, min_dist=0.5, n_neigh=10):
    #pca = PCA(n_components=30)
    #pca.fit(representations)
    #X_reduced = pca.transform(representations)
    #X_reduced = np.multiply(X_reduced, pca.explained_variance_ratio_)

    embedding = umap.UMAP(n_neighbors=n_neigh,
                      min_dist=min_dist,
                      metric='euclidean')
    u = embedding.fit_transform(representations)
    #umap.plot.points(u)
    #u = embedding.fit_transform(X_reduced)
    
    #------------------------------------------------------------------------------
    CAND_COLORS = cm.rainbow(np.linspace(0, 1, 2))
    legend_elements = [Line2D([0], [0], marker='o', color='w', label='GEX',
                          markerfacecolor=CAND_COLORS[0], markersize=8),
                        Line2D([0], [0], marker='o', color='w', label='ATAC',
                          markerfacecolor=CAND_COLORS[1], markersize=8)]
    domains = np.zeros(len(u))
    domains[len(u)//2:] = 1
    domains = domains.astype('int')
    plt.clf()
    fig, ax = plt.subplots(figsize=(5,5))
    scatter = ax.scatter(u[:, 0], u[:, 1], s=1, alpha=0.5, c=CAND_COLORS[domains])
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    #legend1 = ax.legend(*scatter.legend_elements(),
    #                loc="upper left", title="Classes")
    #ax.add_artist(legend1)
    ax.legend(handles=legend_elements, loc='upper left')
    plt.show() #savefig('eval/%s/umap_domains_%s_%s.png'%(expname, lamb, margin))
    
    
    #------------------------------------------------------------------------------
    cell_type = np.unique(labels).tolist()
    nlabels = len(cell_type)
    CAND_COLORS = cm.rainbow(np.linspace(0, 1, nlabels))
    legend_elements = []
    for i in range(nlabels):
        legend_elements.append(Line2D([i],[i],marker='o', color='w', label='%s'%(ct_dic[i]), markerfacecolor=CAND_COLORS[i], markersize=6))
        
    plt.clf()
    fig, ax = plt.subplots(figsize=(5,5))
    ax.scatter(u[:, 0], u[:, 1], s=1, alpha=0.5, c=CAND_COLORS[labels])
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1, 1))
    plt.show() #savefig('eval/%s/umap_labels_%s_%s.png'%(expname, lamb, margin))
    
    #------------------------------------------------------------------------------
    batch_list = np.unique(batch).tolist()
    nlabels = len(batch_list)
    CAND_COLORS = cm.rainbow(np.linspace(0, 1, nlabels+2))
    legend_elements = []
    for i in range(nlabels):
        legend_elements.append(Line2D([i],[i],marker='o', color='w', label='%s'%(batch_dic[i]), markerfacecolor=CAND_COLORS[i+1], markersize=6))
    
    plt.clf()
    fig, ax = plt.subplots(figsize=(5,5))
    ax.scatter(u[:, 0], u[:, 1], s=1, alpha=0.5, c=CAND_COLORS[batch+1])
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.legend(handles=legend_elements, loc='upper left')
    plt.show()
    


    
def umap_plot_split(representations, labels, ct_dic, batch, batch_dic, min_dist=0.5, n_neigh=10):
    #pca = PCA(n_components=30)
    #pca.fit(representations)
    #X_reduced = pca.transform(representations)
    #X_reduced = np.multiply(X_reduced, pca.explained_variance_ratio_)

    embedding = umap.UMAP(n_neighbors=n_neigh,
                      min_dist=min_dist,
                      metric='euclidean')
    u = embedding.fit_transform(representations)
    #u = embedding.fit_transform(X_reduced)
    
    #------------------------------------------------------------------------------
    CAND_COLORS = cm.rainbow(np.linspace(0, 1, 2))
    legend_elements = [Line2D([0], [0], marker='o', color='w', label='GEX',
                          markerfacecolor=CAND_COLORS[0], markersize=8),
                        Line2D([0], [0], marker='o', color='w', label='ATAC',
                          markerfacecolor=CAND_COLORS[1], markersize=8)]
    domains = np.zeros(len(u))
    domains[len(u)//2:] = 1
    domains = domains.astype('int')
    plt.clf()
    fig, ax = plt.subplots(figsize=(5,5))
    scatter = ax.scatter(u[0:len(u)//2, 0], u[0:len(u)//2, 1], s=3, alpha=0.5, c=CAND_COLORS[0])
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    #legend1 = ax.legend(*scatter.legend_elements(),
    #                loc="upper left", title="Classes")
    #ax.add_artist(legend1)
    ax.legend(handles=legend_elements, loc='upper left')
    plt.show() #savefig('eval/%s/umap_domains_%s_%s.png'%(expname, lamb, margin))
    
    
    domains = np.zeros(len(u))
    domains[len(u)//2:] = 1
    domains = domains.astype('int')
    plt.clf()
    fig, ax = plt.subplots(figsize=(5,5))
    scatter = ax.scatter(u[len(u)//2:, 0], u[len(u)//2:, 1], s=3, alpha=0.5, c=CAND_COLORS[1])
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    #legend1 = ax.legend(*scatter.legend_elements(),
    #                loc="upper left", title="Classes")
    #ax.add_artist(legend1)
    ax.legend(handles=legend_elements, loc='upper left')
    plt.show()
    
    
    #------------------------------------------------------------------------------
    batch_list = np.unique(batch).tolist()
    nlabels = len(batch_list)
    CAND_COLORS = cm.rainbow(np.linspace(0, 1, nlabels+2))
    legend_elements = []
    for i in range(nlabels):
        legend_elements.append(Line2D([i],[i],marker='o', color='w', label='%s'%(batch_dic[i]), markerfacecolor=CAND_COLORS[i+1], markersize=6))
    
    plt.clf()
    fig, ax = plt.subplots(figsize=(5,5))
    u_1 = u[batch==0]
    ax.scatter(u_1[:, 0], u_1[:, 1], s=3, alpha=0.5, c=CAND_COLORS[1])
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.legend(handles=legend_elements, loc='upper left')
    plt.show()
    
    plt.clf()
    fig, ax = plt.subplots(figsize=(5,5))
    u_2 = u[batch==1]
    ax.scatter(u_2[:, 0], u_2[:, 1], s=3, alpha=0.5, c=CAND_COLORS[2])
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.legend(handles=legend_elements, loc='upper left')
    plt.show()
    
    