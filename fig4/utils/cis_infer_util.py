import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse

import scglue
from scipy.cluster.hierarchy import ward, dendrogram, leaves_list
from scipy.spatial.distance import pdist



def get_gloc_from_atac_data(peaks):
    """
    Method to get the genomic locations (including the middle point)of peaks
    Author: Linhua Wang Linhua.Wang@bcm.edu
    https://github.com/LiuzLab/Neurips2021/blob/master/task1_utils1.py
    """
    glocs = peaks.tolist()
    glocs = [c for c in glocs if 'chr' in c]
    chrms, ranges, sts, ends, midpoints = [], [], [], [], []
    for gl in glocs:
        chrms.append(gl.split('-')[0])
        st, end = int(gl.split('-')[1]), int(gl.split('-')[2])
        sts.append(st)
        ends.append(end)
        midpoints.append(int((st + end)/2))
        ranges.append("_".join(gl.split("-")[1:]))
    gloc_df = pd.DataFrame({'chrm': chrms, 'grange': ranges,
                        'start': sts, 'end': ends,
                        'midpoint': midpoints}, index=glocs)
    return gloc_df




def get_close_peaks_for_gene(gene, gloc_df, ref_tss_fn="hg38_ref_TSS.txt",
                                                up_dist=150000, down_dist=0):
    """
    Method to get peaks that are within certain range of the TSS site of a gene
    Author: Linhua Wang Linhua.Wang@bcm.edu
    https://github.com/LiuzLab/Neurips2021/blob/master/task1_utils1.py
    """
    
    ref_tss = pd.read_csv(ref_tss_fn, sep='\t')
    ents = ref_tss.loc[ref_tss.symbol == gene]

    if ents.empty:
        return pd.DataFrame()

    peak_ids = []
    peak_dists = []
    promoters = []

    for i in range(len(ents)):
        chrm = ents.seqnames.tolist()[i]
        st = ents.start.tolist()[i]
        end = ents.end.tolist()[i]
        strand = ents.strand.tolist()[i]
        chrm_gloc = gloc_df.loc[gloc_df.chrm == chrm, :]
        if strand == '+':
            peaks = chrm_gloc.loc[((chrm_gloc['midpoint'] >= st - up_dist) & (chrm_gloc['midpoint'] <= end)) |
                    ((chrm_gloc['midpoint'] <= end + down_dist) & (chrm_gloc['midpoint'] >= end))]
            dists = (peaks.midpoint - st).tolist()

        else:
            peaks = chrm_gloc.loc[((chrm_gloc['midpoint'] <= end + up_dist) & (chrm_gloc['midpoint'] >= st)) |
                ((chrm_gloc['midpoint'] >= st - down_dist) & (chrm_gloc['midpoint'] <= st))]
            dists = (end  - peaks.midpoint).tolist()
        peak_ids += peaks.index.tolist()
        peak_dists += dists
        promoters += [1 if (d >= -2000 and d <= 0) else 0 for d in dists]
    df = pd.DataFrame({"peaks": peak_ids, "tss_dist": peak_dists, "pRegion": promoters})
    df['genes'] = gene
    df = df.loc[:,['genes','peaks','tss_dist','pRegion']]
    return df





def empirical_pseudo(anndat, n_pseudocells=250):
    cell_types = anndat.obs['cell_type'].unique()
    pseudo_id = np.zeros(int(anndat.n_obs))
    pseudo_num = anndat.n_obs // n_pseudocells
    for i in range(len(cell_types)):
        celltype_n = int(np.sum(anndat.obs['cell_type']==cell_types[i]))
        n_psedo = celltype_n // pseudo_num + 1
        ct_pseudo_id = np.zeros(celltype_n)
        for j in range(n_psedo-1):
            ct_pseudo_id[(pseudo_num*j):(pseudo_num*(j+1))] = j+1

        ct_pseudo_id[((n_psedo-1)*pseudo_num):] = n_psedo
        #ct_pseudo_id = ct_pseudo_id.astype(int).astype(str)
        #print(ct_pseudo_id.shape)
        #print((anndat.obs['cell_type']==cell_types[i]).shape)
        #print(pseudo_id.shape)
        pseudo_id[anndat.obs['cell_type']==cell_types[i]] = ct_pseudo_id

    return pseudo_id



def hierarchical_based_pseudo(anndat, meta_size = 20):
    cell_type_list = anndat.obs['cell_type'].unique().tolist()
    pseudocell_id = np.zeros(anndat.n_obs)
    ids = 1

    for i in range(len(cell_type_list)):
        num_cells_ct = np.sum(anndat.obs['cell_type'] == cell_type_list[i])
        index_celltype = np.arange(anndat.n_obs)[anndat.obs['cell_type'] == cell_type_list[i]]
        n = 0
        while n<num_cells_ct:
            end = n + meta_size
            if end > num_cells_ct:
                end = num_cells_ct
            pseudocell_id[index_celltype[n:end]] = ids
            ids += 1
            n = end


    assert np.min(pseudocell_id) == 1
    return pseudocell_id





def pseudo_cell(rna, atac, use_rep='SiaNN', n_pseudocells=100, X_agg='mean'):
    print("Clustering pseudocells...")
    if use_rep == 'paired':
        print("Paired mode engaged, clustering based on RNA dataset...")
        Z = ward(pdist(rna.obsm['X_umap']))
        cell_order = leaves_list(Z)
        rna = rna[cell_order]
        atac = atac[cell_order]
        meta_size = rna.n_obs // n_pseudocells
        rna.obs["pseudocell"] = hierarchical_based_pseudo(rna, meta_size = meta_size)
        atac.obs["pseudocell"] = rna.obs["pseudocell"].to_numpy()
        
    elif use_rep == 'SiaNN':
        print("Representation mode engaged, clustering based on combined %s embedding..."%(use_rep))
        combined = anndata.AnnData(
            obs=pd.concat([rna.obs.loc[:,['cell_type']], atac.obs.loc[:,['cell_type']]], join="inner"),
            obsm={use_rep: np.concatenate([rna.obsm[use_rep], atac.obsm[use_rep]])}
        )
        Z = ward(pdist(combined.obsm[use_rep]))
        cell_order = leaves_list(Z)
        combined = combined[cell_order]
        meta_size = rna.n_obs // n_pseudocells
        pseudocell = hierarchical_based_pseudo(combined, meta_size = meta_size)

        combined.obs.insert(1,'pseudocell',pseudocell)
        
        rna.obs["pseudocell_id"] = combined.obs.loc[rna.obs_names, "pseudocell"]
        rna.obs["pseudocell"] = rna.obs["pseudocell_id"] #.astype(str) + rna.obs["cell_type"].astype(str)
        atac.obs["pseudocell_id"] = combined.obs.loc[atac.obs_names, "pseudocell"]
        atac.obs["pseudocell"] = atac.obs["pseudocell_id"] #.astype(str) + atac.obs["cell_type"].astype(str)

    elif use_rep == 'random':
        print("Random mode engaged, clustering based on cell types...")
        ids_seq = np.arange(n_pseudocells)
        n_seq = rna.n_obs // n_pseudocells + 1
        ids_seq = np.repeat(ids_seq,n_seq)
        ids_seq = ids_seq[0:rna.n_obs]
        ids_seq = ids_seq[np.random.choice(rna.n_obs, rna.n_obs, replace=False)]
        rna.obs["pseudocell"] = ids_seq
        
        ids_seq = np.arange(n_pseudocells)
        n_seq = atac.n_obs // n_pseudocells + 1
        ids_seq = np.repeat(ids_seq,n_seq)
        ids_seq = ids_seq[0:atac.n_obs]
        ids_seq = ids_seq[np.random.choice(atac.n_obs, atac.n_obs, replace=False)]
        atac.obs["pseudocell"] = ids_seq
        
        
    elif use_rep == 'empirical':
        print("Empirical mode engaged, clustering based on combined dataset...")
        rna.obs["pseudocell_id"] = empirical_pseudo(rna, n_pseudocells)
        rna.obs["pseudocell"] = rna.obs["cell_type"].astype(str) + rna.obs["pseudocell_id"].astype(str)
        #print(rna.obs["pseudocell"].value_counts())
        atac.obs["pseudocell_id"] = empirical_pseudo(atac, n_pseudocells)
        atac.obs["pseudocell"] = atac.obs["cell_type"].astype(str) + atac.obs["pseudocell_id"].astype(str)
        #print(atac.obs["pseudocell"].value_counts())

    pseudo_label_index = rna.obs[['pseudocell','cell_type']]
    rna_agg = scglue.data.aggregate_obs(rna, "pseudocell", X_agg=X_agg)
    atac_agg = scglue.data.aggregate_obs(atac, "pseudocell", X_agg=X_agg)
    common_pseudocells = np.intersect1d(rna_agg.obs_names, atac_agg.obs_names)
    rna_agg = rna_agg[common_pseudocells, :]
    atac_agg = atac_agg[common_pseudocells, :]

    sc.pp.normalize_total(rna_agg)
    sc.pp.log1p(rna_agg)
    sc.pp.normalize_total(atac_agg)
    sc.pp.log1p(atac_agg)

    return rna_agg, atac_agg, pseudo_label_index





def peaks_within_distance(genes, peaks, distance, ref_tss_fn):
    gloc_df = get_gloc_from_atac_data(peaks)
    df = get_close_peaks_for_gene(genes[0], gloc_df, ref_tss_fn=ref_tss_fn, up_dist=distance)
    for gene in genes[1:]:
        df_tmp = get_close_peaks_for_gene(gene, gloc_df, ref_tss_fn=ref_tss_fn, up_dist=distance)
        df = pd.concat([df,df_tmp] ,ignore_index=True)

    return df





def cis_reg_cor(rna_pseudo, atac_pseudo, peaks_nearby, nonzero=True):
    if nonzero:
        print('Zeros are taken as dropout and removed from pearson correlation calculation.')
        rna_pseudo.X[rna_pseudo.X.toarray()==0] = None
        atac_pseudo.X[atac_pseudo.X.toarray()==0] = None

        nonezero_size = pd.DataFrame(
            ((~np.isnan(rna_pseudo.X.toarray())).T * 1) @ (~np.isnan(atac_pseudo.X.toarray()) * 1),
            index=rna_pseudo.var_names, columns=atac_pseudo.var_names
            )

        rna_m = pd.DataFrame(rna_pseudo.X.toarray())
        atac_m = pd.DataFrame(atac_pseudo.X.toarray())
        corm = pd.concat([rna_m,atac_m],axis=1)
        pcc_score_matrix = corm.corr(method='pearson')
        pcc_score_matrix = pcc_score_matrix.iloc[0:rna_m.shape[1],rna_m.shape[1]:]
        pcc_score_matrix.index = rna_pseudo.var_names
        pcc_score_matrix.columns = atac_pseudo.var_names

        #pcc_score_matrix = pd.DataFrame(
        #    scglue.num.pcc_mat(rna_pseudo.X, atac_pseudo.X),
        #    u.corr(v)
        #    index=rna_pseudo.var_names, columns=atac_pseudo.var_names
        #)

        peaks_nearby['Pearson.cor'] = 0
        for gene in peaks_nearby['genes'].tolist():
            peaks_tmp = peaks_nearby.loc[peaks_nearby['genes']==gene,'peaks'].to_numpy()
            cors_tmp = pcc_score_matrix.loc[gene,peaks_tmp].to_numpy()
            peaks_nearby.loc[peaks_nearby['genes']==gene,'Pearson.cor'] = cors_tmp

            size_tmp = nonezero_size.loc[gene,peaks_tmp].to_numpy()
            peaks_nearby.loc[peaks_nearby['genes']==gene,'test.size'] = size_tmp

        return peaks_nearby

    else:
        print('Zeros are kept and spearman correlation calculation will be done.')
        spr_score_matrix = pd.DataFrame(
            scglue.num.spr_mat(rna_pseudo.X, atac_pseudo.X),
            index=rna_pseudo.var_names, columns=atac_pseudo.var_names
        )


        #print(spr_score_matrix)
        peaks_nearby['Spearman.cor'] = 0
        for gene in peaks_nearby['genes'].tolist():
            peaks_tmp = peaks_nearby.loc[peaks_nearby['genes']==gene,'peaks'].to_numpy()
            cors_tmp = spr_score_matrix.loc[gene,peaks_tmp].to_numpy()
            #print(peaks_nearby.loc[peaks_nearby['genes']==gene,'Spearman.cor'])
            #print(cors_tmp)
            peaks_nearby.loc[peaks_nearby['genes']==gene,'Spearman.cor'] = cors_tmp

        return peaks_nearby







def cis_element_score(rna, atac, genes, distance, 
                      ref_tss_fn, use_rep, n_pseudocells=200,
                      X_agg='mean', nonzero=False,
                      return_pseudo_bulk=False):
    """
    inferring gene's cis-regulatory element from its nearby peaks
    Parameters
    ----------
    rna
        snRNA-seq AnnData dataset
    atac
        snATAC-seq AnnData dataset
    genes
        list of genes to be inferred for cis-regulatory element
    distance
        distances from gene TSS to be considered

    Returns
    -------
    cis_score_df
        cis-regulatory element score for listed genes
    """
    print('Add suffix to cell names in RNA and ATAC data...')
    rna.obs_names = rna.obs_names + '_RNA'
    atac.obs_names = atac.obs_names + '_ATAC'

    rna_pseudo, atac_pseudo, pseudo_label_index = pseudo_cell(rna, atac, use_rep, n_pseudocells=n_pseudocells, X_agg=X_agg)
    #print(rna_pseudo)
    #print(atac_pseudo)
    rna_pseudo = rna_pseudo[:, genes]

    peaks = atac_pseudo.var_names
    peaks_nearby = peaks_within_distance(genes, peaks, distance, ref_tss_fn)
    #print(peaks_nearby)
    
    print('Selecting peaks within %sbp of the genes...'%(distance))
    peaks = peaks_nearby['peaks'].unique()
    atac_pseudo = atac_pseudo[:, peaks]

    cis_score_df = cis_reg_cor(rna_pseudo, atac_pseudo, peaks_nearby, nonzero=nonzero)
    
    if return_pseudo_bulk:
        return cis_score_df, rna_pseudo, atac_pseudo, pseudo_label_index
    else:
        return cis_score_df


