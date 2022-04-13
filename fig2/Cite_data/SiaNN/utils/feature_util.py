import anndata as ad
import numpy as np
import scipy as sp
import sys
import os


#meta = {'resources_dir': '.'}
#sys.path.append(f"{meta['resources_dir']}")


def feature_selection_cite(input_mod1, input_mod2, path):
    gene_index = open('%s/gene_index_cite.txt'%(path),'r').readlines()
    gene_index = [i.strip() for i in gene_index]
    gene_index = np.array(gene_index)
    prt_index = open('%s/protein_index.txt'%(path),'r').readlines()
    prt_index = [i.strip() for i in prt_index]
    prt_index = np.array(prt_index)

    intersect_gene = np.intersect1d(input_mod1.var_names.to_numpy(), gene_index)
    cpl_gene = np.setdiff1d(gene_index, input_mod1.var_names.to_numpy())
    if len(cpl_gene) > 0:
        print('Complementing missing feaures with 0')
        input_mod1 = input_mod1[:,intersect_gene]
        cpl_matrix = sp.sparse.csr_matrix((input_mod1.n_obs, len(cpl_gene)), dtype=np.float64)
        input_mod1_cpl = ad.AnnData(
            X=cpl_matrix,
            uns={
                "dataset_id": input_mod1.uns['dataset_id'],
            }
        )

        input_mod1_cpl.var_names = cpl_gene
        input_mod1_cpl.obs_names = input_mod1.obs_names

        test = ad.concat([input_mod1, input_mod1_cpl], axis=1)
        test.obs = input_mod1.obs
        test.uns = input_mod1.uns
        test = test[:,gene_index]
        input_mod1 = test
    else:
        input_mod1 = input_mod1[:,gene_index]


    intersect_prt = np.intersect1d(input_mod2.var_names.to_numpy(), prt_index)
    cpl_prt = np.setdiff1d(prt_index, input_mod2.var_names.to_numpy())
    if len(cpl_prt) > 0:
        print('Complementing missing feaures with 0')
        input_mod2 = input_mod2[:,intersect_prt]
        cpl_matrix = sp.sparse.csr_matrix((input_mod2.n_obs, len(cpl_prt)), dtype=np.float64)
        input_mod2_cpl = ad.AnnData(
            X=cpl_matrix,
            uns={
                "dataset_id": input_mod2.uns['dataset_id'],
            }
        )

        input_mod2_cpl.var_names = cpl_prt
        input_mod2_cpl.obs_names = input_mod2.obs_names

        test = ad.concat([input_mod2, input_mod2_cpl], axis=1)
        test.obs = input_mod2.obs
        test.uns = input_mod2.uns
        test = test[:,prt_index]
        input_mod2 = test
    else:
        input_mod2 = input_mod2[:,prt_index]

    return input_mod1, input_mod2










