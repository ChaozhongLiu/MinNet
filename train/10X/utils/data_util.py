from tqdm import tqdm, trange
import statistics
import math
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sp
import scanpy as sc

from matplotlib import colors
from sklearn.preprocessing import normalize

import os
import anndata
import logging

from sklearn.preprocessing import normalize



class dataset(object):
    def __init__(self, mod1, mod2, batch_number=250, validation=True):
        super(dataset, self).__init__()
        self.mod1 = mod1
        self.mod2 = mod2
        self.batch_number = batch_number
        self.validation = validation
        #self.anchor_score = anchor_score
        self.init_dataset()



    def init_dataset(self):
        # load dataset and save as self.data_set_X and self.data_set_Y
        print('Loading mod1 from %s'%self.mod1)
        self.input_mod1 = ad.read_h5ad(self.mod1)

        print('Loading mod2 from %s'%self.mod2)
        self.input_mod2 = ad.read_h5ad(self.mod2)

        print('mod1: ', self.input_mod1.X.shape)
        print('mod2: ', self.input_mod2.X.shape)

        # making cell order the same
        self.input_mod2 = self.input_mod2[self.input_mod1.obs_names]
        assert np.all(self.input_mod1.obs_names == self.input_mod1.obs_names)

        # make gene and peak order the same
        print('Selecting and ordering features...')
        gene_index = open('data/gene_index.txt','r').readlines()
        gene_index = [i.strip() for i in gene_index]
        gene_index = np.array(gene_index)
        self.input_mod1 = self.input_mod1[:,gene_index]
        self.input_mod2 = self.input_mod2[:,gene_index]

        # normalizing RNA-seq data
        print('Normalizing data...')
        sc.pp.scale(self.input_mod1, max_value=10)
        sc.pp.scale(self.input_mod2, max_value=10)
        self.input_mod1.X = sp.sparse.csr_matrix(self.input_mod1.X)
        self.input_mod2.X = sp.sparse.csr_matrix(self.input_mod2.X)


        # batch information
        self.input_mod1.obs['batch'] = 'sd'
        self.input_mod2.obs['batch'] = 'sd'
        self.batch_list = self.input_mod1.obs['batch'].unique().tolist()
        self.batch_num = len(self.batch_list) # not batch_number for training size
        print('%s batches detected, will train the model with batch separated in each epoch.'%(self.batch_num))

        # cell type information
        self.cell_type_name = self.input_mod1.obs['cell_type'].unique().to_list()
        self.input_mod1.obs['cell_type_name'] = self.input_mod1.obs['cell_type']
        self.input_mod2.obs['cell_type_name'] = self.input_mod2.obs['cell_type']
        self.input_mod1.obs['cell_type'] = self.input_mod1.obs['cell_type'].cat.codes
        self.input_mod2.obs['cell_type'] = self.input_mod1.obs['cell_type']
        self.label_set = self.input_mod1.obs['cell_type'].unique().tolist()
        self.label_num = len(self.label_set)
        print('%s labels detected.'%(self.label_num))
        self.label_size = self.input_mod1.obs.groupby(["batch", "cell_type"]).size().reset_index(name="size")

        # label modality
        self.input_mod1.obs['acc'] = 0
        self.input_mod2.obs['acc'] = 1



        #split into train and validation
        if self.validation:
            print('Spliting each batch into training and validation...')
            self.n_cells = len(self.input_mod1.obs_names)
            self.num_val = self.n_cells // 5
            self.num_train = self.n_cells - self.num_val

            # split train and test
            np.random.seed(2021)
            index_train = np.random.choice(self.n_cells, self.num_train, replace=False)
            self.index_train = np.sort(index_train, axis=None)
            mask = np.zeros(len(self.input_mod1.obs_names), np.bool)
            mask[index_train] = 1
            self.input_mod1.obs['train'] = mask
            self.input_mod2.obs['train'] = mask
            self.input_mod1_val = self.input_mod1[self.input_mod1.obs['train']==0]
            self.input_mod2_val = self.input_mod2[self.input_mod2.obs['train']==0]
            self.input_mod1 = self.input_mod1[self.input_mod1.obs['train']==1]
            self.input_mod2 = self.input_mod2[self.input_mod2.obs['train']==1]

            assert np.all(self.input_mod1.obs_names == self.input_mod2.obs_names)
            assert np.all(self.input_mod1_val.obs_names == self.input_mod2_val.obs_names)


            #batch size
            self.batch_size = self.num_train//self.batch_number
            self.batch_size_val = self.num_val//self.batch_number

            print('Training dataset:')
            print('\t%s cells,'%(int(len(self.input_mod1.obs_names))))

            print()
            print('Validation dataset:')
            print('\t%s cells,'%(int(len(self.input_mod1_val.obs_names))))
            

        else:
            print('No validation set')
            self.n_cells = len(self.input_mod1.obs_names)
            self.num_cells = len(self.input_mod1.obs_names)
            self.num_train = len(self.input_mod1.obs_names)

            self.input_mod1.obs['train'] = 1
            self.input_mod2.obs['train'] = 1

            self.batch_size = self.num_train//self.batch_number

            print('Training dataset:')
            print('\t%s cells,'%(int(len(self.input_mod1.obs_names))))
            print()




    def prepare_siamese(self, mode='default'):

        self.input_mod1_pair = -np.ones(self.num_train, dtype=int)
        self.input_mod2_pair = -np.ones(self.num_train, dtype=int)

        self._train_z_mod1 = -np.ones(self.num_train, dtype=int)
        self._train_z_mod2 = -np.ones(self.num_train, dtype=int)

        dis_index_mod1 = -np.ones(self.num_train, dtype=int)
        dis_index_mod2 = -np.ones(self.num_train, dtype=int)


        for index in range(self.num_train):
            negative_ids = np.where(np.arange(self.num_train) != index)[0]

            # prepare mod1 pair in mod2
            if np.random.random() > 0.75:
                self.input_mod1_pair[index] = index
                self._train_z_mod1[index] = 1
                dis_index_mod1[index] = index
                #self._use_domain_mod1[index] = 0

            elif len(negative_ids) >= 1:
                siamese_index = np.random.choice(negative_ids)
                self.input_mod1_pair[index] = siamese_index
                self._train_z_mod1[index] = 0
                dis_index_mod1[index] = siamese_index
                #self._use_domain_mod1[index] = self.input_mod1.obsp['distance'][index,siamese_index]
            #else:
            #    self._use_domain_mod1[index]  = 0

            # prepare mod2 pair in mod1
            if np.random.random() > 0.75:
                self.input_mod2_pair[index] = index
                self._train_z_mod2[index] = 1
                dis_index_mod2[index] = index
                #self._use_domain_mod2[index] = 0
            elif len(negative_ids) >= 1:
                siamese_index = np.random.choice(negative_ids)
                self.input_mod2_pair[index] = siamese_index
                self._train_z_mod2[index] = 0
                dis_index_mod2[index] = siamese_index
                #self._use_domain_mod2[index] = self.input_mod1.obsp['distance'][index,siamese_index]

        self._use_domain_mod1 = self.input_mod1.obsp['distance'][np.arange(self.num_train),dis_index_mod1]
        self._use_domain_mod2 = self.input_mod1.obsp['distance'][np.arange(self.num_train),dis_index_mod2]
        
        ct_diff_index_mod1 = 3.0 * (self.input_mod1.obs['cell_type'].to_numpy() != self.input_mod1.obs['cell_type'][dis_index_mod1].to_numpy())
        ct_diff_index_mod2 = 3.0 * (self.input_mod1.obs['cell_type'].to_numpy() != self.input_mod1.obs['cell_type'][dis_index_mod2].to_numpy())
        
        self._use_domain_mod1 += ct_diff_index_mod1
        self._use_domain_mod2 += ct_diff_index_mod2


        assert np.all(self.input_mod1_pair != -1)
        assert np.all(self.input_mod2_pair != -1)
        assert np.all(self._train_z_mod1 != -1)
        assert np.all(self._train_z_mod2 != -1)

        if self.validation:

            self.input_mod1_pair_val = -np.ones(self.num_val, dtype=int)
            self.input_mod2_pair_val = -np.ones(self.num_val, dtype=int)

            self._val_z_mod1 = -np.ones(self.num_val, dtype=int)
            self._val_z_mod2 = -np.ones(self.num_val, dtype=int)

            #self._use_domain_mod1_val = np.ones(self.num_val)
           # self._use_domain_mod2_val = np.ones(self.num_val)
            dis_index_mod1 = -np.ones(self.num_val, dtype=int)
            dis_index_mod2 = -np.ones(self.num_val, dtype=int)

            for index in range(self.num_val):
                negative_ids = np.where(np.arange(self.num_val) != index)[0]

                # prepare mod1 pair in mod2
                if np.random.random() > 0.75:
                    self.input_mod1_pair_val[index] = index
                    self._val_z_mod1[index] = 1
                    dis_index_mod1[index] = index
                    #self._use_domain_mod1_val[index] = 0
                elif len(negative_ids) >= 1:
                    siamese_index = np.random.choice(negative_ids)
                    self.input_mod1_pair_val[index] = siamese_index
                    self._val_z_mod1[index] = 0
                    dis_index_mod1[index] = siamese_index
                    #self._use_domain_mod1_val[index] = self.input_mod1_val.obsp['distance'][index,siamese_index]


                # prepare mod1 pair in mod2
                if np.random.random() > 0.75:
                    self.input_mod2_pair_val[index] = index
                    self._val_z_mod2[index] = 1
                    dis_index_mod2[index] = index
                    #self._use_domain_mod2_val[index] = 0
                elif len(negative_ids) >= 1:
                    siamese_index = np.random.choice(negative_ids)
                    self.input_mod2_pair_val[index] = siamese_index
                    self._val_z_mod2[index] = 0
                    dis_index_mod2[index] = siamese_index
                    #self._use_domain_mod2_val[index] = self.input_mod1_val.obsp['distance'][index,siamese_index]

            self._use_domain_mod1_val = self.input_mod1_val.obsp['distance'][np.arange(self.num_val),dis_index_mod1]
            self._use_domain_mod2_val = self.input_mod1_val.obsp['distance'][np.arange(self.num_val),dis_index_mod2]
            
            ct_diff_index_mod1 = 3.0 * (self.input_mod1_val.obs['cell_type'].to_numpy() != self.input_mod1_val.obs['cell_type'][dis_index_mod1].to_numpy())
            ct_diff_index_mod2 = 3.0 * (self.input_mod1_val.obs['cell_type'].to_numpy() != self.input_mod1_val.obs['cell_type'][dis_index_mod2].to_numpy())

            self._use_domain_mod1_val += ct_diff_index_mod1
            self._use_domain_mod2_val += ct_diff_index_mod2



    
    def train_data(self, loss=0.0, intervals=2, epoch=0, reorder=True):

        if (epoch % int(intervals) == 0) & (epoch % 10 != 4):
            print('Generating siamese pairs...')
            self.prepare_siamese(mode='default')
        elif (epoch % int(intervals) == 0) & (epoch % 10 == 4):
            #print('Generating siamese pairs with highly similar negative pairs...')
            #self.prepare_siamese(mode='hard')
            print('Generating siamese pairs...')
            self.prepare_siamese(mode='default')


        if reorder:
            rand_index = np.random.choice(self.num_train, self.num_train, replace=False)
            self.rand_index = rand_index
        else:
            rand_index = np.arange(self.num_train)
            self.rand_index = rand_index
        #self.input_mod1_tmp = self.input_mod1_tmp[rand_index]
        self.input_mod1_pair_tmp = self.input_mod1_pair[rand_index]
        self._train_z_mod1_tmp = self._train_z_mod1[rand_index]
        self._use_domain_mod1_tmp = self._use_domain_mod1[rand_index]

        #self.input_mod2_tmp = self.input_mod2_tmp[rand_index]
        self.input_mod2_pair_tmp = self.input_mod2_pair[rand_index]
        self._train_z_mod2_tmp = self._train_z_mod2[rand_index]
        self._use_domain_mod2_tmp = self._use_domain_mod2[rand_index]


        return tqdm(self.next_train_batch(),
                    desc='Train loss: {:.4f}'.format(loss),
                    total=self.batch_number, mininterval=1.0, leave=False)



    
    def next_train_batch(self):
        start = 0
        #startY = 0
        end = self.batch_size #originally end
        #endY = self.batch_size_Y #add new
        N = self.num_train
        #NY = len(self._train_Y_tmp)

        while start < N:

            if self.validation:
                start_val = np.random.randint(0, self.num_val-self.batch_size)
                end_val = start_val + self.batch_size
                val_mod1 = self.input_mod1_val[start_val:end_val].X.toarray()
                val_mod1_pair = self.input_mod2_val[self.input_mod1_pair_val[start_val:end_val]].X.toarray()
                val_label = self.input_mod1_val[start_val:end_val].obs['cell_type'].to_numpy()
                val_mod2 = self.input_mod2_val[start_val:end_val].X.toarray()
                val_mod2_pair = self.input_mod1_val[self.input_mod2_pair_val[start_val:end_val]].X.toarray()


            if end < N:
                train_mod1 = self.input_mod1[self.rand_index[start:end]].X.toarray()
                train_mod2 = self.input_mod2[self.rand_index[start:end]].X.toarray()
                train_mod1_pair = self.input_mod2[self.input_mod1_pair_tmp[start:end]].X.toarray()
                train_mod2_pair = self.input_mod1[self.input_mod2_pair_tmp[start:end]].X.toarray()
                train_label = self.input_mod1[self.rand_index[start:end]].obs['cell_type'].to_numpy()

                if self.validation:
                    yield train_mod1, train_mod2, train_mod1_pair, train_mod2_pair, train_label,\
                        self._train_z_mod1_tmp[start:end], self._train_z_mod2_tmp[start:end],\
                        self._use_domain_mod1_tmp[start:end], self._use_domain_mod2_tmp[start:end],\
                        val_mod1, val_mod2, val_mod1_pair, val_mod2_pair, val_label,\
                        self._val_z_mod1[start_val:end_val], self._val_z_mod2[start_val:end_val],\
                        self._use_domain_mod1_val[start_val:end_val], self._use_domain_mod2_val[start_val:end_val]
                else:
                    yield train_mod1, train_mod2, train_mod1_pair, train_mod2_pair, train_label,\
                        self._train_z_mod1_tmp[start:end], self._train_z_mod2_tmp[start:end],\
                        self._use_domain_mod1_tmp[start:end], self._use_domain_mod2_tmp[start:end],\
                        None, None, None, None, None, None, None, None, None
            
            else:
                train_mod1 = self.input_mod1[self.rand_index[start:]].X.toarray()
                train_mod2 = self.input_mod2[self.rand_index[start:]].X.toarray()
                train_mod1_pair = self.input_mod2[self.input_mod1_pair_tmp[start:]].X.toarray()
                train_mod2_pair = self.input_mod1[self.input_mod2_pair_tmp[start:]].X.toarray()
                train_label = self.input_mod1[self.rand_index[start:]].obs['cell_type'].to_numpy()

                if self.validation:
                    yield train_mod1, train_mod2, train_mod1_pair, train_mod2_pair, train_label,\
                        self._train_z_mod1_tmp[start:], self._train_z_mod2_tmp[start:],\
                        self._use_domain_mod1_tmp[start:], self._use_domain_mod2_tmp[start:],\
                        val_mod1, val_mod2, val_mod1_pair, val_mod2_pair, val_label,\
                        self._val_z_mod1[start_val:end_val], self._val_z_mod2[start_val:end_val],\
                        self._use_domain_mod1_val[start_val:end_val], self._use_domain_mod2_val[start_val:end_val]
                else:
                    yield train_mod1, train_mod2, train_mod1_pair, train_mod2_pair, train_label,\
                        self._train_z_mod1_tmp[start:], self._train_z_mod2_tmp[start:],\
                        self._use_domain_mod1_tmp[start:], self._use_domain_mod2_tmp[start:],\
                        None, None, None, None, None, None, None, None, None

            start = end
            end += self.batch_size






class test_data(object):
    def __init__(self, mod1, mod2, batch_number=250):
        super(test_data, self).__init__()
        self.mod1 = mod1
        self.mod2 = mod2
        self.batch_number = batch_number
        self.init_dataset()


    def init_dataset(self):
        # load dataset and save as self.data_set_X and self.data_set_Y
        print('Loading mod1 from %s'%self.mod1)
        self.input_mod1 = ad.read_h5ad(self.mod1)

        print('Loading mod2 from %s'%self.mod2)
        self.input_mod2 = ad.read_h5ad(self.mod2)

        print('mod1: ', self.input_mod1.X.shape)
        print('mod2: ', self.input_mod2.X.shape)

        # making cell order the same
        self.input_mod2 = self.input_mod2[self.input_mod1.obs_names]
        assert np.all(self.input_mod1.obs_names == self.input_mod1.obs_names)

        # make gene and peak order the same
        print('Selecting and ordering features...')
        gene_index = open('data/gene_index.txt','r').readlines()
        gene_index = [i.strip() for i in gene_index]
        gene_index = np.array(gene_index)
        self.input_mod1 = self.input_mod1[:,gene_index]
        self.input_mod2 = self.input_mod2[:,gene_index]


        # normalizing RNA-seq data with 0 mean and 1 sd
        print('Normalizing data...')
        sc.pp.scale(self.input_mod1, max_value=10)
        sc.pp.scale(self.input_mod2, max_value=10)
        self.input_mod1.X = sp.sparse.csr_matrix(self.input_mod1.X)
        self.input_mod2.X = sp.sparse.csr_matrix(self.input_mod2.X)


        # label modality
        self.input_mod1.obs['acc'] = 0
        self.input_mod2.obs['acc'] = 1

        #self.n_cells = len(self.input_mod1.obs_names)
        self.num_cells = len(self.input_mod1.obs_names)
        self.batch_size = self.num_cells//self.batch_number


        


