import torch
import torch.nn as nn
import numpy as np
import math
from torch.autograd import Variable



class SiamAtt(nn.Module):

    def __init__(self, num_peaks, num_genes, label_num, dim1=768, dim2=32):
        super().__init__()
        self.num_peaks = num_peaks
        self.num_genes = num_genes
        
        self.dim1 = dim1
        self.dim2 = dim2
        
        #----------- Mod 2 -------------
        # Creat input layer
        self.MLP_mod2 = nn.Sequential(
            nn.Linear(self.num_genes, self.dim1),
            nn.Tanh(),
            nn.BatchNorm1d(self.dim1),
        )
        

        #----------- Mod 1 -------------
        # Creat MLP to project num_channels to num_genes
        self.MLP_mod1 = nn.Sequential(
            nn.Linear(self.num_genes, self.dim1),
            nn.Tanh(),
            nn.BatchNorm1d(self.dim1),
        )
        
        
        #----------- Loss module -------------
        # Desciminator of the Siamese network
        self.domain_mlp = nn.Sequential(
            nn.Linear(self.dim1, self.dim2),
        )

        self.classifier = nn.Sequential(
            nn.Linear(self.dim1, label_num),
        )

        #print(self)


        
    def forward(self, x, y, mode='train'):

        y = self.MLP_mod2(y)
             
        x = self.MLP_mod1(x)

        class_x = self.classifier(x)
        class_y = self.classifier(y)

        domain_x = self.domain_mlp(x)
        domain_y = self.domain_mlp(y)
        if mode == 'train':
            return domain_x, domain_y, class_x, class_y
        elif mode == 'test':
            return x, y, class_x, class_y




