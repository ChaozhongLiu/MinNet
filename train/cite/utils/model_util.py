import torch
import torch.nn as nn
import numpy as np
import math
from torch.autograd import Variable


class SiamAtt(nn.Module):

    def __init__(self, num_prts, num_genes, label_num, dim1=768, dim2=32):
        super().__init__()
        self.num_prts = num_prts
        self.num_genes = num_genes
        
        self.dim1 = dim1
        self.dim2 = dim2


        # -------------------------------------------------------------------- #
        #  processing of GEX data
        # -------------------------------------------------------------------- #
        # Creat input layer
        self.MLP_mod1 = nn.Sequential(
            nn.Linear(self.num_genes, self.dim1),
            nn.Sigmoid(),
            nn.BatchNorm1d(self.dim1),
        )
        

        # -------------------------------------------------------------------- #
        #  processing of ADT data
        # -------------------------------------------------------------------- #
        # Creat MLP to project num_channels to num_genes
        self.MLP_mod2 = nn.Sequential(
            nn.Linear(self.num_prts, self.dim1),
            nn.Sigmoid(),
            nn.BatchNorm1d(self.dim1),
        )
        

        # -------------------------------------------------------------------- #
        #  loss modules
        # -------------------------------------------------------------------- #
        self.domain_mlp = nn.Sequential(
            nn.Linear(self.dim1, self.dim2),
        )

        self.classifier = nn.Sequential(
            nn.Linear(self.dim1, label_num),
        )

        #print(self)


        
    def forward(self, x, y, mode='train'):

        #--------- gex ----------
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

        return domain_x, domain_y, class_x, class_y







