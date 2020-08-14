# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import cPickle as pickle

from collections import OrderedDict

pd.set_option('io.hdf.default_format','table')

# Global path
path_to_target_db = '../data/targets/'

# Local paths
paths = OrderedDict()
paths['coexpression'] = 'coexpression/GeneExp_AllHuman.txt'
paths['essentiality'] = 'essentiality/Essentiality_AllHuman.txt'
paths['functional_association'] = 'functional_association/FxnAss_AllHuman.txt'
paths['genetic_interaction'] = 'genetic_interactions/GeneticInteraction_AllHuman.txt'
paths['go_all_shared_fraction'] = 'go_all_shared_fraction/Fraction_All.txt'
paths['go_bp_shared_fraction'] = 'go_bp_shared_fraction/Fraction_BP.txt'
paths['go_cc_shared_fraction'] = 'go_cc_shared_fraction/Fraction_CC.txt'
paths['go_mf_shared_fraction'] = 'go_mf_shared_fraction/Fraction_MF.txt'
paths['go_slim_shared_fraction'] = 'go_slim_shared_fraction/Fraction_Slim.txt'
paths['pathway_shared_fraction'] = 'pathway_shared_fraction/Pathway_Fraction_All.txt'
paths['shortest_path_length'] = 'shortest_path_length/ShortestPathLength_AllHuman.txt'


# Either make a list of all human protein names from the gene expression data
if True:
    set_of_protein_names = set()
    with open(path_to_target_db + paths['coexpression'], 'r') as fh:
        for line_number, line in enumerate(fh):
            try:
                row = [l.strip() for l in line.split('\t')]
                protein_1, protein_2 = row[1].split('_')
                set_of_protein_names.add(protein_1)
                set_of_protein_names.add(protein_2)
                if line_number % 10000 == 0:
                    print line_number
            except ValueError:
                print line
                raise
else:
    with open('../data/set_of_protein_names.pickle', 'r') as fh:
        set_of_protein_names = pickle.load(fh)
    list_of_protein_names = sorted(list(set_of_protein_names))
        
"""
data_shape = (len(list_of_protein_names),len(list_of_protein_names))

def symmetrize(a):
    return a + a.T - np.diag(a.diagonal())

def fill_data(data_type):
    data = pd.DataFrame()
    if data_type == 'essentiality' or data_type == 'genetic_interactions':
        protein_column_idx = 0
    else:
        protein_column_idx = 1
    with open(path_to_target_db + paths[data_type]) as fh:
        for line_number, line in enumerate(fh):
            if line_number % 1000 == 0:
                print line_number
                print line
            break
            row = [l.strip() for l in line.split('\t')]
            protein_1, protein_2 = sorted(row[protein_column_idx].split('_'))
            value = float(row[2])
            data.ix[protein_1, protein_2] = value
#    data = symmetrize(data)
    return pd.DataFrame(data=data, index=list_of_protein_names, columns=list_of_protein_names)
            

store = pd.HDFStore('../data/data.h5', 'w')
for data_type in paths:
    temp_df = fill_data(data_type)
    print temp_df
    print temp_df.shape
    store[data_type] = temp_df
#data_type = paths.keys()[0]

#temp_df = fill_data(data_type)
#store[data_type] = temp_df


#temp_df_1 = temp_df[:5]
#temp_df_2 = temp_df[8:12]



#store.append('coexpression',temp_df_1)
#store.append('coexpression',temp_df_2)

store.close()

"""