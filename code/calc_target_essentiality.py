# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd

class essentiality_db():
    def __init__(self, input_essentiality_file):
        
        with open(input_essentiality_file, 'rb') as fh:
            data = [ [row.split('\t')[0].strip(), row.split('\t')[1].strip()] for row in fh ]
            
        self.target_essentiality = dict( (a,b) for a, b in data )


    def get_pair(self, ensp_pair):
        ensp1 = ensp_pair.split('_')[0].strip()
        ensp2 = ensp_pair.split('_')[1].strip()
        
        if not (self.target_essentiality.has_key(ensp1) and self.target_essentiality.has_key(ensp2)):
            return 0
        
        try:
            gene_ess1 = int(self.target_essentiality[ensp1])
            gene_ess2 = int(self.target_essentiality[ensp2])
        except ValueError:
            return 0
            
        feature = [gene_ess1, gene_ess2]
        if feature == [-1,-1]:
            result = 1
        elif feature == [-1,0] or feature == [0,-1]:
            result = 2
        elif feature == [0,0]:
            result = 3
        elif feature == [1,0] or feature == [0,1]:
            result = 4
        elif feature == [1,1]:
            result = 5
        elif feature == [-1,1] or feature == [1,-1]:
            result = 6
        else:
            raise Exception("Didn't take something into account!")
        return result

if __name__ == '__main__':
    from class_db import project_db

    data_path = '../data/'
    
    p_db = project_db()
    
    coexpr_db = essentiality_db(data_path + 'targets/Essentiality/Human_GeneEss_ENSP.txt')    
    
    target_pairs_pos = np.array(p_db.target_data_pos.index)
    essentiality_data_pos = np.vectorize(coexpr_db.get_pair)(target_pairs_pos)
    essentiality_pos_df = pd.DataFrame(data=essentiality_data_pos, index=target_pairs_pos, columns=['essentiality'])
    with open(data_path + 'targets/Essentiality/essentiality_pos.tsv', 'w') as fh:
        essentiality_pos_df.to_csv(fh, sep='\t')
    
    
    target_pairs_neg = np.array(p_db.target_data_neg.index)
    essentiality_data_neg = np.vectorize(coexpr_db.get_pair)(target_pairs_neg)
    essentiality_neg_df = pd.DataFrame(data=essentiality_data_neg, index=target_pairs_neg, columns=['essentiality'])
    with open(data_path + 'targets/Essentiality/essentiality_neg.tsv', 'w') as fh:
        essentiality_neg_df.to_csv(fh, sep='\t')
    