# -*- coding: utf-8 -*-

"""
based on the file by fei
"""

import numpy as np
import pandas as pd
from collections import defaultdict, OrderedDict


class ATCDB():

    def __init__(self, atc_dict):
        self.atc_dict = atc_dict


    def _split_atc_into_levels(self, atc_code):
        return atc_code[0], atc_code[1:3], atc_code[3], atc_code[4], atc_code[5:7]


    def get_sim(self, cid_1, cid_2):
        if not (self.atc_dict.has_key(cid_1) and self.atc_dict.has_key(cid_2)):
            print 'No atc data found for drugs with cids: {}, {}'.format(cid_1, cid_2)
            return np.nan

        highest_score = 0
        for code1 in self.atc_dict[cid_1]:
            split_code1 = self._split_atc_into_levels(code1)

            for code2 in self.atc_dict[cid_2]:
                split_code2 = self._split_atc_into_levels(code2)

                score = 0
                if (split_code1[0] == split_code2[0]):
                    score = 1
                    if (split_code1[1] == split_code2[1]):
                        score = 2
                        if (split_code1[2] == split_code2[2]):
                            score = 3
                            if (split_code1[3] == split_code2[3]):
                                score = 4
                                if (split_code1[4] == split_code2[4]):
                                    score = 5
                if (score > highest_score):
                   highest_score = score

        return highest_score



class SideEffectDB():

    def __init__(self, side_effect_dict):
        self.side_effect_dict = side_effect_dict


    def get_sim(self, cid_1, cid_2):
        if not (self.side_effect_dict.has_key(cid_1) and self.side_effect_dict.has_key(cid_2)):
            print 'No side effect data found for drugs with cids: {}, {}'.format(cid_1, cid_2)
            return np.nan

        set_1 = self.side_effect_dict[cid_1]
        set_2 = self.side_effect_dict[cid_2]
        score = len(set_1 & set_2)/float(len(set_1 | set_2))
        return score




if __name__ == '__main__':
    path_to_data = '/home/kimlab1/strokach/working/databases/chemical_interactions/version_2.1/'

    def calculate_similarity(filename, DB, data_column_name):
        with open(path_to_data + filename) as fh:
            raw_data = fh.readlines()
        data_dict = defaultdict(set)
        for line in raw_data:
            cid, values = line.strip().split('\t')
            cid = int(cid)
            values = values.replace('[','').replace(']','').replace("'",'').replace(' ','').split(',')
            data_dict[cid].update(values)
        db = DB(data_dict)
        data_dict_keys = data_dict.keys()
        data_dict_keys.sort()
        final_data = OrderedDict([('cid_1', []), ('cid_2', []), (data_column_name, [])])
        for idx_1 in range(len(data_dict_keys)):
            for idx_2 in range(idx_1, len(data_dict_keys)):
                cid_1 = data_dict_keys[idx_1]
                cid_2 = data_dict_keys[idx_2]
                final_data['cid_1'].append(cid_1)
                final_data['cid_2'].append(cid_2)
                final_data[data_column_name].append(db.get_sim(cid_1, cid_2))
        final_data_df = pd.DataFrame(final_data)
        final_data_df.to_csv(path_to_data + filename + '.with_sim_score', sep='\t', index=False)

    calculate_similarity('drug_atc.tsv', ATCDB, 'atc_similarity')
    calculate_similarity('drug_side_effect.tsv', SideEffectDB, 'side_effect_similarity')


# Obsolete:
#    from class_db import project_db
#
#    data_path = '../data/'
#
#    p_db = project_db()
#    my_atc_db = atc_db(p_db.drug_ATC)
#
#    drug_pairs_pos = np.array(p_db.drug_data_pos.index)
#    atc_data_pos = np.vectorize(my_atc_db.get_sim)(drug_pairs_pos)
#    atc_data_pos_df = pd.DataFrame(data=atc_data_pos, index=drug_pairs_pos, columns=['atc_similarity'])
#    with open(data_path + 'drugs/clinical-similarity/atc_similarity_pos.tsv', 'w') as fh:
#        atc_data_pos_df.to_csv(fh, sep='\t')
#
#    drug_pairs_neg = np.array(p_db.drug_data_neg.index)
#    atc_data_neg = np.vectorize(my_atc_db.get_sim)(drug_pairs_neg)
#    atc_data_neg_df = pd.DataFrame(data=atc_data_neg, index=drug_pairs_neg, columns=['atc_similarity'])
#    with open(data_path + 'drugs/clinical-similarity/atc_similarity_neg.tsv', 'w') as fh:
#        atc_data_neg_df.to_csv(fh, sep='\t')