import torch
import torch.nn as nn
import numpy as np
import math
from torch.autograd import Variable



class SiamAtt_multiome(nn.Module):

    def __init__(self, num_peaks, num_genes, channel_start=16, channels=256, num_conv=2, 
                 num_att=4, num_heads=8, dropout=0.2, domain_dim1=32, domain_dim2=32):
        super().__init__()
        self.num_peaks = num_peaks
        self.num_genes = num_genes
        #self.channel_start = channel_start
        self.channels = channels
        #self.num_conv = num_conv
        #self.num_att = num_att
        #self.num_heads = num_heads
        #self.dropout = dropout
        self.domain_dim1 = domain_dim1
        #self.domain_dim2 = domain_dim2
        #self.num_chunks = self.num_peaks
        
        #----------- Mod 2 -------------
        # Creat input layer
        self.MLP_mod2 = nn.Sequential(
            nn.Linear(self.num_genes, self.channels*3),
            nn.Tanh(),
            nn.BatchNorm1d(self.channels*3),
        )


        #----------- Mod 1 -------------
        # Creat MLP to project num_channels to num_genes
        self.MLP_mod1 = nn.Sequential(
            nn.Linear(self.num_genes, self.channels*3),
            nn.Tanh(),
            nn.BatchNorm1d(self.channels*3),
        )
        
        
        #----------- Loss module -------------
        # Desciminator of the Siamese network
        #self.mlp_bn_mod1 = nn.BatchNorm1d(filter_list_2[-1])
        #self.mlp_bn_mod2 = nn.BatchNorm1d(filter_list_2[-1])
        self.domain_mlp = nn.Sequential(
            #nn.BatchNorm1d(self.channels),
            nn.Linear(self.channels*3, self.domain_dim1),
            #nn.Tanh(),
            #nn.Linear(self.domain_dim1, self.domain_dim2),
            #nn.Tanh()
        )

        self.classifier = nn.Sequential(
            nn.Linear(self.channels*3, 22),
            #nn.Linear(filter_list_2[-1], 21),
            #nn.Sigmoid(),
        )

        #print(self)


        
    def forward(self, x, y, mode='matching'):

        #input
        y = self.MLP_mod2(y)
        
        x = self.MLP_mod1(x)

        class_x = self.classifier(x)
        class_y = self.classifier(y)

        domain_x = self.domain_mlp(x)
        domain_y = self.domain_mlp(y)

        if mode == 'matching':
            return domain_x, domain_y, class_x, class_y
        elif mode == 'embedding':
            return x, y, class_x, class_y






