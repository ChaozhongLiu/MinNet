import torch
import torch.nn as nn
import numpy as np
import math
from torch.autograd import Variable




#--------------
class SiamAtt_cite(nn.Module):

    def __init__(self, num_prts, num_genes, channel_start=64, channels=256, num_conv=4, 
                 num_att=4, num_heads=8, dropout=0.2, domain_dim1=32):
        super().__init__()
        self.num_prts = num_prts
        self.num_genes = num_genes
        #self.channel_start = channel_start
        self.channels = channels
        #self.num_conv = num_conv
        #self.num_att = num_att
        #self.num_heads = num_heads
        #self.dropout = dropout
        self.domain_dim1 = domain_dim1
        #self.domain_dim2 = domain_dim2
        #self.num_chunks = self.num_genes


        # -------------------------------------------------------------------- #
        #  processing of GEX data
        # -------------------------------------------------------------------- #
        # Creat input layer
        self.MLP_mod1 = nn.Sequential(
            #nn.LayerNorm(self.channels),
            #nn.BatchNorm1d(self.num_genes),
            nn.Linear(self.num_genes, self.channels*3),
            nn.Sigmoid(),
            #nn.Linear(self.channels*6, self.channels*3),
            #nn.Sigmoid(),
            nn.BatchNorm1d(self.channels*3),
            #nn.Linear(filter_list_2[-1]*2, self.domain_dim1),
            #nn.Tanh(),
        )


        # -------------------------------------------------------------------- #
        #  processing of ADT data
        # -------------------------------------------------------------------- #
        # Creat MLP to project num_channels to num_genes
        self.MLP_mod2 = nn.Sequential(
            #nn.LayerNorm(self.channels),
            #nn.BatchNorm1d(self.num_genes),
            nn.Linear(self.num_prts, self.channels*3),
            nn.Sigmoid(),
            #nn.Linear(self.channels*6, self.channels*3),
            #nn.Sigmoid(),
            nn.BatchNorm1d(self.channels*3),
            #nn.Linear(filter_list_2[-1]*2, self.domain_dim1),
            #nn.Tanh(),
        )
        


        # -------------------------------------------------------------------- #
        #  loss modules
        # -------------------------------------------------------------------- #
        self.domain_mlp = nn.Sequential(
            #nn.BatchNorm1d(self.channels),
            nn.Linear(self.channels*3, self.domain_dim1),
            #nn.Tanh(),
            #nn.Linear(self.domain_dim1, self.domain_dim2),
            #nn.Tanh()
        )

        self.classifier = nn.Sequential(
            nn.Linear(self.channels*3, 42),
            #nn.Linear(filter_list_2[-1], 21),
            #nn.Sigmoid(),
        )

        #print(self)


        
    def forward(self, x, y, mode='matching'):

        #--------- gex ----------
        #input
        x = self.MLP_mod1(x)

        #--------- adt ----------
        y = self.MLP_mod2(y)

        #--------- loss ----------
        class_x = self.classifier(x)
        class_y = self.classifier(y)
        #pred = self.mlp_bn_1(pred)
        #x = self.mlp_bn_2(x)
        domain_x = self.domain_mlp(x)
        domain_y = self.domain_mlp(y)

        if mode == 'matching':
            return domain_x, domain_y, class_x, class_y
        elif mode == 'embedding':
            return x, y, class_x, class_y



