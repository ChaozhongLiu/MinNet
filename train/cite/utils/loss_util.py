import torch.nn as nn
import torch.nn.functional as F

class ContrastiveLoss(nn.Module):
    """
    Contrastive loss
    Takes embeddings of two samples and a target label == 1 if samples are from the same class and label == 0 otherwise
    """
    def __init__(self, margin):
        super(ContrastiveLoss, self).__init__()
        self.margin = margin
    def forward(self, output1, output2, target, use_domain, size_average=True):
        #distances = F.cosine_embedding_loss(output1, output2, target, margin=self.margin, reduction='none')
        distances = (output2 - output1).pow(2).sum(1)+0.001  # squared distances
        losses = 0.5 * (target.float() * distances +
                        2 * (1 + -1 * target).float() * F.relu(self.margin*use_domain - distances.sqrt()).pow(2))
        #losses = losses*use_domain
        #print(distances[target < 1].sqrt())
        return losses.mean() if size_average else losses
