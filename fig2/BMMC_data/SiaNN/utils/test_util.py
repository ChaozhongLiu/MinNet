import anndata as ad
import numpy as np
import pandas as pd
import scipy as sp
import scanpy as sc

import torch
import torch.nn as nn
import torch.optim as optim
from torch.autograd import Variable
import sys


from utils.model_util_test import *


# To load the dataset
class test_data_multiome(object):
    def __init__(self, mod1, mod2, batch_number=200):
        super(test_data_multiome, self).__init__()
        self.input_mod1 = mod1
        self.input_mod2 = mod2
        self.batch_number = batch_number
        self.init_dataset()


    def init_dataset(self):

        # label modality
        self.input_mod1.obs['acc'] = 0
        self.input_mod2.obs['acc'] = 1

        #self.n_cells = len(self.input_mod1.obs_names)
        self.num_cells_1 = len(self.input_mod1.obs_names)
        self.num_cells_2 = len(self.input_mod2.obs_names)
        self.batch_size_1 = self.num_cells_1//self.batch_number
        self.batch_size_2 = self.num_cells_2//self.batch_number





# To setup the network for test
class Siamese_Test_multiome(object):
    def __init__(self, dataset, num_peaks, num_genes):
        super(Siamese_Test_multiome, self).__init__()
        self.dataset = dataset
        self.D = SiamAtt_multiome(num_peaks, num_genes).cuda()


    def initialize(self):
        for p in self.D.parameters():
            if len(p.shape) > 1:
                nn.init.xavier_uniform(p)
            else:
                nn.init.normal(p, 0, 1)

    def load_model(self, model_file):
        skpt_dict = torch.load(model_path)
        self.D.load_state_dict(skpt_dict)

    def save_model(self, model_file):
        torch.save(self.D.state_dict(), model_file)



    def test(self):
        #self.D.train()
        self.D.eval()

        output_X = None
        output_Y = None

        n_iter = self.dataset.num_cells_1 // self.dataset.batch_size_1
        rng_state = np.random.get_state()
        for i in range(n_iter):
            x = self.dataset.input_mod1[i*self.dataset.batch_size_1:(i+1)*self.dataset.batch_size_1].X.toarray()
            y = self.dataset.input_mod2[i*self.dataset.batch_size_2:(i+1)*self.dataset.batch_size_2].X.toarray()

            X = Variable(torch.cuda.FloatTensor(x))
            Y = Variable(torch.cuda.FloatTensor(y))

            output_x, output_y, class_x, class_y = self.D(X, Y)

            if output_X is None:
                output_X = output_x.cpu().data.numpy()
                output_Y = output_y.cpu().data.numpy()
                class_X = class_x.cpu().data.numpy()
                class_Y = class_y.cpu().data.numpy()
            else:
                output_X = np.concatenate((output_X, output_x.cpu().data.numpy()), 0)
                output_Y = np.concatenate((output_Y, output_y.cpu().data.numpy()), 0)
                class_X = np.concatenate((class_X, class_x.cpu().data.numpy()), 0)
                class_Y = np.concatenate((class_Y, class_y.cpu().data.numpy()), 0)
            del X, Y, output_x, output_y

        # last batch
        x = self.dataset.input_mod1[(i+1)*self.dataset.batch_size_1:].X.toarray()
        y = self.dataset.input_mod2[(i+1)*self.dataset.batch_size_2:].X.toarray()

        X = Variable(torch.cuda.FloatTensor(x))
        Y = Variable(torch.cuda.FloatTensor(y))

        output_x, output_y, class_x, class_y = self.D(X, Y)

        output_X = np.concatenate((output_X, output_x.cpu().data.numpy()), 0)
        output_Y = np.concatenate((output_Y, output_y.cpu().data.numpy()), 0)
        class_X = np.concatenate((class_X, class_x.cpu().data.numpy()), 0)
        class_Y = np.concatenate((class_Y, class_y.cpu().data.numpy()), 0)

        return output_X, output_Y, class_X, class_Y


