from torch.autograd import Function
import torch.nn.functional as F
import pandas as pd
import numpy as np
import torch
import matplotlib.pyplot as plt
from torch.autograd import Variable
import matplotlib.cm as cm
from matplotlib.patches import Patch
from matplotlib.lines import Line2D



def split_loss(output1, output2, target, use_domain):

    distances = (output2 - output1).pow(2).sum(1)+0.001
    loss_pos = ((distances*use_domain)[(target == 1)]).mean()
    loss_neg = ((distances*use_domain)[(target == 0)]).mean()
    loss_pos_list = distances[(target == 1)]
    loss_neg_list = distances[(target == 0)]

    return loss_pos, loss_neg, loss_pos_list, loss_neg_list



def loss_hist(pos_loss, neg_loss, pos_loss_val, neg_loss_val):
    kwargs = dict(alpha=0.5, bins=100, density=False, stacked=True)
    plt.hist(pos_loss, **kwargs, color='dodgerblue', label='train_pos')
    plt.hist(neg_loss, **kwargs, color='orange', label='train_neg')
    plt.gca().set(title='Histogram of Distances between Siamese Pairs - Train', xlabel='Distance')
    #plt.xlim(50,75)
    plt.legend()
    plt.show()
    plt.hist(pos_loss_val, **kwargs, color='deeppink', label='val_pos')
    plt.hist(neg_loss_val, **kwargs, color='green', label='val_neg')
    plt.gca().set(title='Histogram of Distances between Siamese Pairs - Val', xlabel='Distance')
    #plt.xlim(50,75)
    plt.legend()
    plt.show()



