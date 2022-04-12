import os
import time
import math
import numpy as np
from tqdm import tqdm

import torch
import torch.nn as nn
import torch.optim as optim
from torch.autograd import Variable

from .model_util import *
from .loss_util import ContrastiveLoss
from .function_util import get_distance, loss_hist, split_loss
import matplotlib.pyplot as plt
from sklearn.metrics import pairwise_distances

from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import normalize


class BaseTrainer(object):
    def __init__(self, num_epochs, batch_number, model_path, use_gpu, validation=True):
        super(BaseTrainer, self).__init__()
        self.num_epochs = num_epochs
        self.batch_number = batch_number
        self.use_gpu = use_gpu
        self.validation = validation
        self.out_path = model_path
        
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

    '''
    def valid(self):
        self.D.eval()
        counts = 0
        sum_acc = 0.0
        for x, y in self.dataset.valid_data():
            counts += len(y)
            if self.use_gpu:
                X = Variable(torch.cuda.FloatTensor(x))
                Y = Variable(torch.cuda.LongTensor(y))
            else:
                X = Variable(torch.FloatTensor(x))
                Y = Variable(torch.LongTensor(y))
            f_X = self.D(X, mode='test')
            sum_acc += (f_X.max(dim=1)[1] == Y).float().sum().data.cpu().numpy()
            del X,Y,f_X
        valid_acc = sum_acc/counts
        print('Valid Accuracy: {0:0.3f}'.format(valid_acc))
        return valid_acc
    '''


class Siamese_Trainer(BaseTrainer):
    def __init__(self, dataset, num_epochs, batch_number, model_path, num_prts, num_genes,
                 margin=10.0, lamb=0.5, use_gpu=True, validation=True, lr=1e-3):
        #d_dim: number of genes; margin:margin of contrastive loss; 
        #Setup network
        super(Siamese_Trainer, self).__init__(num_epochs, batch_number, model_path, use_gpu, validation)
        self.dataset = dataset
        self.margin = margin
        self.lamb = lamb
        self.lr = lr
        self.D = SiamAtt(num_prts, num_genes, int(self.dataset.label_num)).cuda()
        self.L = ContrastiveLoss(self.margin).cuda()
        self.C = nn.CrossEntropyLoss().cuda()

        self.optimizer = optim.Adam([{'params':self.D.parameters()}], lr=self.lr, weight_decay=1e-6)
        
        
    def train(self, f, initial = True, plots=True, intervals=2):
        self.D.train()
        if initial:
            self.initialize()
        siam_loss_val = float('inf')
        self.train_acc_list = []
        self.train_losses = [] #training classification loss
        #self.adv_loss = [] #training adversarial loss
        self.loss_pos = []
        self.loss_neg = []
        self.valid_accs = [] #validation accuracy
        self.valid_losses = [] #validation adversarial loss
        best_validate_loss = float('inf')


        for j in range(self.num_epochs):
            begin = time.time()
            self.D.train()
            
            counts = 0
            #counts_x = 0
            counts_val = 0
            sum_acc = 0.0 
            val_sum_acc = 0.0
            train_epoch_loss = []
            #adv_epoch_loss = []
            epoch_loss_pos = []
            epoch_loss_neg = []
            valid_epoch_loss = []
            valid_epoch_domain = []
            pos_list = []
            neg_list = []
            #val_list = []
            pos_list_val = []
            neg_list_val = []
            train_data = self.dataset.train_data(float(siam_loss_val), intervals=intervals, epoch = j)
            for x1, y1, x2, y2, label, zx, zy, dx, dy, x1val, y1val, x2val, y2val, label_val, zxval, zyval, dxval, dyval in train_data:
                #forward calculation and back propagation
                X1 = Variable(torch.cuda.FloatTensor(x1))
                Y1 = Variable(torch.cuda.FloatTensor(y1))
                X2 = Variable(torch.cuda.FloatTensor(x2))
                Y2 = Variable(torch.cuda.FloatTensor(y2))
                X = torch.cat((X1,Y2), dim=0)
                Y = torch.cat((X2,Y1), dim=0)
                Label = Variable(torch.cuda.LongTensor(label))
                #Label = torch.cat((Label,Label), dim=0)
                ZX = Variable(torch.cuda.FloatTensor(zx))
                ZY = Variable(torch.cuda.FloatTensor(zy))
                Z = torch.cat((ZX,ZY), dim=0)
                DX = Variable(torch.cuda.FloatTensor(dx))
                DY = Variable(torch.cuda.FloatTensor(dy))
                D = torch.cat((DX,DY), dim=0)

                self.optimizer.zero_grad()
                output_x, output_y, class_x, class_y = self.D(X,Y)
                domain_loss = self.L(output_x, output_y, Z, D)

                num_c = len(class_x) // 2
                class_x = class_x[0:num_c]
                class_y = class_y[num_c:]
                label_loss = self.C(class_x, Label) + self.C(class_y, Label)

                # Save the loss for plot
                pos_loss, neg_loss, loss_pos_list, loss_neg_list = split_loss(output_x, output_y, Z, D)
                pos_loss = pos_loss.data.cpu().numpy()
                neg_loss = neg_loss.data.cpu().numpy()
                loss_pos_list = loss_pos_list.data.cpu().numpy().tolist()
                pos_list += loss_pos_list
                loss_neg_list = loss_neg_list.data.cpu().numpy().tolist()
                neg_list += loss_neg_list

                # Regulation
                #penalty_loss_x = torch.cat([x.view(-1) for x in self.D.feature_extractor_X.parameters()])
                #penalty_loss_y = torch.cat([x.view(-1) for x in self.D.feature_extractor_Y.parameters()])
                #penalty_loss_x = lambda1 * torch.norm(penalty_loss_x, 1)
                #penalty_loss_y = lambda2 * torch.norm(penalty_loss_y, 1)

                self.loss = label_loss + self.lamb * domain_loss
                siam_loss_val = self.loss.data.cpu().numpy()
                train_data.set_description('Train loss: {0:.4f}'.format(float(siam_loss_val)))
                self.loss.backward()
                self.optimizer.step()

                train_epoch_loss.append(siam_loss_val)
                #adv_epoch_loss.append(domain_loss_val)
                epoch_loss_pos.append(pos_loss)
                epoch_loss_neg.append(neg_loss)

                #record accuracy
                counts += x1.shape[0] + y1.shape[0]
                sum_acc += (class_x.max(dim=1)[1] == Label).float().sum().data.cpu().numpy() + (class_y.max(dim=1)[1] == Label).float().sum().data.cpu().numpy()

                del X1, Y1, X2, Y2, Label, X, Y, ZX, ZY, Z, DX, DY, D, output_x, output_y, class_x, class_y, label_loss, domain_loss, pos_loss, neg_loss, loss_pos_list, loss_neg_list

                if self.validation:
                    X1v = Variable(torch.cuda.FloatTensor(x1val))
                    Y1v = Variable(torch.cuda.FloatTensor(y1val))
                    X2v = Variable(torch.cuda.FloatTensor(x2val))
                    Y2v = Variable(torch.cuda.FloatTensor(y2val))
                    X_v = torch.cat((X1v,Y2v), dim=0)
                    Y_v = torch.cat((X2v,Y1v), dim=0)
                    Label_v = Variable(torch.cuda.LongTensor(label_val))
                    #Label = torch.cat((Label,Label), dim=0)
                    ZX = Variable(torch.cuda.FloatTensor(zxval))
                    ZY = Variable(torch.cuda.FloatTensor(zyval))
                    Z = torch.cat((ZX,ZY), dim=0)
                    DX = Variable(torch.cuda.FloatTensor(dxval))
                    DY = Variable(torch.cuda.FloatTensor(dyval))
                    D = torch.cat((DX,DY), dim=0)
                    
                    output_x, output_y, class_x, class_y = self.D(X_v,Y_v)
                    domain_loss = self.L(output_x, output_y, Z, D)

                    num_c = len(class_x) // 2
                    class_x = class_x[0:num_c]
                    class_y = class_y[num_c:]
                    label_loss = self.C(class_x, Label_v) + self.C(class_y, Label_v)

                    loss_val = label_loss + self.lamb * domain_loss
                    loss_val = loss_val.data.cpu().numpy()
                    valid_epoch_loss.append(loss_val)

                    _, _, loss_pos_list, loss_neg_list = split_loss(output_x, output_y, Z, D)
                    loss_pos_list = loss_pos_list.data.cpu().numpy().tolist()
                    pos_list_val += loss_pos_list
                    loss_neg_list = loss_neg_list.data.cpu().numpy().tolist()
                    neg_list_val += loss_neg_list

                    #record accuracy
                    counts_val += x1val.shape[0] + y1val.shape[0]
                    val_sum_acc += (class_x.max(dim=1)[1] == Label_v).float().sum().data.cpu().numpy() + (class_y.max(dim=1)[1] == Label_v).float().sum().data.cpu().numpy()

                    del X1v,Y1v,X2v,Y2v,X_v,Y_v,Label_v,ZX,ZY,Z,DX,DY,D,output_x,output_y,class_x,class_y,domain_loss,label_loss,loss_val,loss_pos_list,loss_neg_list
                
            #print(sum_acc, counts)
            self.train_acc = sum_acc/counts #counts
            self.train_acc_list.append(self.train_acc)
            self.train_losses.append(np.sum(train_epoch_loss)/self.batch_number)
            #self.adv_loss.append(np.sum(adv_epoch_loss)/self.batch_number)
            self.loss_pos.append(np.mean(epoch_loss_pos))
            self.loss_neg.append(np.mean(epoch_loss_neg))
            if self.validation:
                valid_acc = val_sum_acc/counts_val
                self.valid_accs.append(valid_acc)
                valid_loss = np.sum(valid_epoch_loss)/self.batch_number
                self.valid_losses.append(valid_loss)

                f.write("Epoch %d, time = %ds, train accuracy = %.4f, train loss = %.4f, validation accuracy = %.4f, validation loss = %.4f\n" % (
                        j, time.time() - begin, self.train_acc, self.train_losses[-1], valid_acc, self.valid_losses[-1]))

                if self.valid_losses[-1] < best_validate_loss:
                    best_validate_loss = self.valid_losses[-1]
                    #score = self.test()
                    f.write("Best Validated Model Loss = %.4f\n" % (best_validate_loss))
                    self.save_model(os.path.join(self.out_path, 'best_model_margin%s_lamb%s.ckpt'%(self.margin, self.lamb)))
            else:
                f.write("Epoch %d, time = %ds, train accuracy = %.4f, train loss = %.4f\n" % (
                    j, time.time() - begin, self.train_acc, self.train_losses[-1]))

            #if j%10==9:
            tqdm.write('After epoch {0} Train Accuracy: {1:0.3f}, Train Loss: {2:0.3f}'.format(j+1, self.train_acc, self.train_losses[-1]))
            if self.validation:
                tqdm.write('After epoch {0} Val Accuracy: {1:0.3f}, Val Loss: {2:0.3f}'.format(j+1, valid_acc, self.valid_losses[-1]))
            
                #output_x, output_y = self.val_test()
                #self.D.train()
                #score_val = self.match_score(output_x, output_y)
                #tqdm.write('After epoch {0} Val Score: {1:0.3f}'.format(j+1, score_val))

            #keep this part currently
            if j+1 in [5,10,15,20,25,30,35,40,45, self.num_epochs] + [i for i in range(50, self.num_epochs, 5)]:

                #if j+1 in [i for i in range(50, self.num_epochs, 25)]:
                self.save_model(os.path.join(self.out_path, 'model_margin%s_lamb%s_epoch%s.ckpt'%(self.margin, self.lamb, j+1)))

                if plots:
                    loss_hist(np.array(pos_list), np.array(neg_list), np.array(pos_list_val), np.array(neg_list_val))

            
        self.save_model(os.path.join(self.out_path, 'final_model_margin%s_lamb%s.ckpt'%(self.margin, self.lamb)))

        #plot learning curve
        if plots:
            plt.figure(figsize=(5,5))
            plt.plot(self.train_losses, label='train_loss')
            if self.validation:
                plt.plot(self.valid_losses, label='val_loss')
            plt.title('Learning Curve - Loss')
            plt.xlabel('epochs')
            plt.legend()
            plt.show()

            plt.figure(figsize=(5,5))
            plt.plot(self.train_acc_list, label='train_acc')
            if self.validation:
                plt.plot(self.valid_accs, label='val_acc')
            plt.title('Learning Curve - Accuracy')
            plt.xlabel('epochs')
            plt.legend()
            plt.show()



    def test(self):
        #self.D.train()
        self.D.eval()

        output_X = None
        output_Y = None

        n_iter = self.dataset.num_cells // self.dataset.batch_size
        rng_state = np.random.get_state()
        for i in range(n_iter):
            x = self.dataset.input_mod1[i*self.dataset.batch_size:(i+1)*self.dataset.batch_size].X.toarray()
            y = self.dataset.input_mod2[i*self.dataset.batch_size:(i+1)*self.dataset.batch_size].X.toarray()

            X = Variable(torch.cuda.FloatTensor(x))
            Y = Variable(torch.cuda.FloatTensor(y))

            output_x, output_y, class_x, class_y = self.D(X, Y, mode='train')

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
            del X, Y, output_x, output_y, class_x, class_y

        # last batch
        x = self.dataset.input_mod1[(i+1)*self.dataset.batch_size:].X.toarray()
        y = self.dataset.input_mod2[(i+1)*self.dataset.batch_size:].X.toarray()

        X = Variable(torch.cuda.FloatTensor(x))
        Y = Variable(torch.cuda.FloatTensor(y))

        output_x, output_y, class_x, class_y = self.D(X, Y, mode='train')

        output_X = np.concatenate((output_X, output_x.cpu().data.numpy()), 0)
        output_Y = np.concatenate((output_Y, output_y.cpu().data.numpy()), 0)
        class_X = np.concatenate((class_X, class_x.cpu().data.numpy()), 0)
        class_Y = np.concatenate((class_Y, class_y.cpu().data.numpy()), 0)

        return output_X, output_Y, class_X, class_Y


    def val_test(self):
        #self.D.train()
        self.D.eval()

        output_X = None
        output_Y = None

        n_iter = self.dataset.num_val // self.dataset.batch_size_val
        rng_state = np.random.get_state()
        for i in range(n_iter):
            x = self.dataset.input_mod1_val[i*self.dataset.batch_size_val:(i+1)*self.dataset.batch_size_val].X.toarray()
            y = self.dataset.input_mod2_val[i*self.dataset.batch_size_val:(i+1)*self.dataset.batch_size_val].X.toarray()

            X = Variable(torch.cuda.FloatTensor(x))
            Y = Variable(torch.cuda.FloatTensor(y))

            output_x, output_y, _, _ = self.D(X, Y)

            if output_X is None:
                output_X = output_x.cpu().data.numpy()
                output_Y = output_y.cpu().data.numpy()
            else:
                output_X = np.concatenate((output_X, output_x.cpu().data.numpy()), 0)
                output_Y = np.concatenate((output_Y, output_y.cpu().data.numpy()), 0)
            del X, Y, output_x, output_y

        # last batch
        x = self.dataset.input_mod1_val[(i+1)*self.dataset.batch_size_val:].X.toarray()
        y = self.dataset.input_mod2_val[(i+1)*self.dataset.batch_size_val:].X.toarray()

        X = Variable(torch.cuda.FloatTensor(x))
        Y = Variable(torch.cuda.FloatTensor(y))

        output_x, output_y, _, _ = self.D(X, Y)

        output_X = np.concatenate((output_X, output_x.cpu().data.numpy()), 0)
        output_Y = np.concatenate((output_Y, output_y.cpu().data.numpy()), 0)

        return output_X, output_Y
