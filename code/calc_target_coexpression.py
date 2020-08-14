# -*- coding: utf-8 -*-

import numpy as np
from scipy.stats import pearsonr
import pandas as pd

class coexpression_db():
    def __init__(self, input_correlation_file):
        
        with open(input_correlation_file, 'rb') as fh:
            coexpression_df = pd.read_csv(fh, sep='\t')
    
        coexpression_index = [l.strip() for l in coexpression_df['SampleCode']]
        del coexpression_df['SampleCode']
        
        coexpression_array = np.array(coexpression_df, dtype=float)
        
        num_proteins = len(coexpression_index)
        pearson_correlation = np.zeros((num_proteins, num_proteins),  float)

        for i in range(0,num_proteins):
            for j in range(0,num_proteins):
                pearson_correlation[i,j] = pearsonr(coexpression_array[i,:],coexpression_array[j,:])[0]
        self.coexpression_index = coexpression_index
        self.pearson_correlation = pearson_correlation


    def get_pearsonr(self,cid):
        print cid
        cid1, cid2 = [c.strip() for c in cid.split('_')]
        print cid1,cid2     
        
        if not (cid1 in self.coexpression_index 
            and cid2 in self.coexpression_index):
            return 0.
        
        bins = [-0.19066943, -0.07752145,  0.01560775,  0.11676259,  0.25336123]
        idx1 = self.coexpression_index.index(cid1)
        idx2 = self.coexpression_index.index(cid2)
        pearson_r = self.pearson_correlation[idx1, idx2]        

        if pearson_r < bins[0]:
            return 1.
        elif pearson_r < bins[1]:
            return 2.
        elif pearson_r < bins[2]:
            return 3.
        elif pearson_r < bins[3]:
            return 4.
        elif pearson_r < bins[4]:
            return 5.
        else:
            return 6.


if __name__ == '__main__':
    from class_db import project_db

    data_path = '../data/'
    
    p_db = project_db()
    
    coexpr_db = coexpression_db(data_path + 'targets/Coexpression/E-MTAB-62_GeneExpTargetProbeNR.txt')    
    
    target_pairs_pos = np.array(p_db.target_data_pos.index)
    coexpression_data_pos = np.vectorize(coexpr_db.get_pearsonr)(target_pairs_pos)
    coexpression_pos_df = pd.DataFrame(data=coexpression_data_pos, index=target_pairs_pos, columns=['coexpression'])
    with open(data_path + 'targets/Coexpression/coexpression_pos.tsv', 'w') as fh:
        coexpression_pos_df.to_csv(fh, sep='\t')
    
    
    target_pairs_neg = np.array(p_db.target_data_neg.index)
    coexpression_data_neg = np.vectorize(coexpr_db.get_pearsonr)(target_pairs_neg)
    coexpression_neg_df = pd.DataFrame(data=coexpression_data_neg, index=target_pairs_neg, columns=['coexpression'])
    with open(data_path + 'targets/Coexpression/coexpression_neg.tsv', 'w') as fh:
        coexpression_neg_df.to_csv(fh, sep='\t')
