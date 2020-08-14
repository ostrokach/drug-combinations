# -*- coding: utf-8 -*-
"""
@author: alexey
"""

import MySQLdb as mdb
import cPickle as pickle
import numpy as np
import pandas as pd
import random
import os
import urllib2 
import json

class project_db():
    
    def __init__(self,
                     drug_pairs_filename='../data/Drug_Pairs_PositiveSet.txt', 
                     drugs_targets_pos_filename='../data/Target_Pairs_Rel2_Positive.txt', 
                     drugs_targets_neg_filename='../data/Target_Pairs_Rel2_Negative.txt'):
        """ 
        Initialise objects to store all the required data
        
        self.drug_pairs (dataframe): 
            'drug_pairs_id': a unique id for each drug pair in the positive set
            'drug_pair': the pubchem (or drugbank for proteins) accession id for the two drugs
        self.drugs_targets_pos (dataframe):
            'drug_pair': the pubchem (or drugbank for proteins) accession id for the two drugs
            'target_pair': the ensembl protein accession for the two drug targets
            # each drug pair should have ~ 4 protein target pairs
        self.drug_targets_neg (dataframe):
            # same as above but for the negative training set (obtained using random pairs of proteins)
        self.drug_data_pos (dataframe):
            Contains data for each feature of the drug-drug pair
            .index
        self.drug_data_neg (dataframe):
            # same as above but for the negative training set
        self.target_data_pos (dataframe):
            Contains data for each feature of the target-target pair
            .index
        self.targets_data_neg (dataframe):
            # same as above but for the negative training set
        self.drug_SMILES (dict):
            Contains the SMILES code for each drug in the positive and negative training set (excluding protein-based drugs)
        """
        
        # DataFrame to convert unique ids Clare assigned to positive drug pairs, to STITCH drug accessions
        with open(drug_pairs_filename, 'r') as ih:
            self.drug_pairs = pd.read_csv(ih, sep='\t', skiprows=1, names=['drug_pair_id', 'drug_pair'], dtype=str)
            self.drug_pairs.set_index('drug_pair_id', drop=False, inplace=True, verify_integrity=False)
        
        # DataFrame to link all drug_pairs with all target_pairs in the positive set
        drugs_targets_pos = set()
        with open(drugs_targets_pos_filename, 'r') as ih:
            for line in ih:
                row = [ l.strip() for l in line.split('\t') ]
                target_pair = row[0]
                drug_pair_ids = row[2].replace('[','').replace(']','').replace("'", '').split(', ')
                for drug_pair_id in drug_pair_ids:
                    drug_pair = self.drug_pairs.loc[drug_pair_id]['drug_pair']
                    if (drug_pair, target_pair,) in drugs_targets_pos:
                        print (drug_pair, target_pair,), 'already in positive set!'
                        continue
                    drugs_targets_pos.add((drug_pair, target_pair,))
        self.drugs_targets_pos = pd.DataFrame(list(drugs_targets_pos),  columns=['drug_pair', 'target_pair'])
        
        # DataFrame to link all drug_pairs with all target_pairs in the negative set
        drugs_targets_neg = set()
        with open(drugs_targets_neg_filename, 'r') as ih:
            for line in ih:
                row = [ l.strip() for l in line.split('\t') ]
                target_pair = row[0]
                drug_pairs = row[2].replace('[','').replace(']','').replace("'", '').split(', ')
                for drug_pair in drug_pairs:
                    if (drug_pair, target_pair,) in drugs_targets_pos:
                        print (drug_pair, target_pair,), 'are also in positive set!'
                        continue
                    if (drug_pair, target_pair,) in drugs_targets_neg:
                        print (drug_pair, target_pair,), 'already in negative set!'
                        continue
                    drugs_targets_neg.add((drug_pair, target_pair,))
        self.drugs_targets_neg = pd.DataFrame(list(drugs_targets_neg), columns=['drug_pair', 'target_pair'])            

 
        #######################################################################       
        # DataFrames to keep all the features for machine learning
        self.drug_data_pos = pd.DataFrame(index=list(set(self.drugs_targets_pos['drug_pair'])))
        self.target_data_pos = pd.DataFrame(index=list(set(self.drugs_targets_pos['target_pair'])))
        self.drug_data_neg = pd.DataFrame(index=list(set(self.drugs_targets_neg['drug_pair'])))
        self.target_data_neg = pd.DataFrame(index=list(set(self.drugs_targets_neg['target_pair'])))


        #######################################################################
        # Get SMILES terms for the positive and negative drug combinations
        list_of_stitch_ids = []
        map(list_of_stitch_ids.extend, [l.split('_') for l in self.drug_data_pos.index])
        map(list_of_stitch_ids.extend, [l.split('_') for l in self.drug_data_neg.index])
        print "Number of redundant drugs (should be > 10000):", len(list_of_stitch_ids)
        list_of_stitch_ids = list(set(list_of_stitch_ids))
        print  "Number of non-redundant drugs:", len(list_of_stitch_ids)
        
        list_of_pubchem_ids = []
        for stitch_id in list_of_stitch_ids:
            if stitch_id[:2] == 'CI':
                list_of_pubchem_ids.append(str(int(stitch_id[4:])))
            elif stitch_id[:2] == 'DB':
                list_of_pubchem_ids.append(stitch_id)

        
        #######################################################################
        # Get the names of compounds       
        if (os.path.isfile('../data/pickled/drug_names.pickle') and
                os.path.isfile('../data/pickled/drug_SMILES.pickle')):
            with open('../data/pickled/drug_names.pickle', 'rb') as fh:
                self.drug_names = pickle.load(fh)
            with open('../data/pickled/drug_SMILES.pickle', 'rb') as fh:
                self.drug_SMILES = pickle.load(fh)
        else:           
            mysql_query = "SELECT chemical, name, SMILES_string \
                            FROM stitch.chemicals \
                            WHERE chemical IN ('%s');" \
                            % "', '".join(list_of_stitch_ids)
            
            con = mdb.connect('localhost', 'root', 'kim630', 'stitch')
            with con:
                cur = con.cursor()
                cur.execute(mysql_query)
                rows = cur.fetchall()
            
            self.drug_names = dict( (a, b) for a, b, c in rows )
            with open('../data/pickled/drug_names.pickle', 'wb') as fh:
                pickle.dump(self.drug_names, fh)         
            self.drug_SMILES = dict( (a, c) for a, b, c in rows )
            with open('../data/pickled/drug_SMILES.pickle', 'wb') as fh:
                pickle.dump(self.drug_SMILES, fh)

        
        #######################################################################
        # Get ATC code for the STITCH database
        if os.path.isfile('../data/pickled/drug_ATC1.pickle'):
            with open('../data/pickled/drug_ATC.pickle', 'rb') as fh:
                self.drug_ATC = pickle.load(fh)
        else:
            mysql_query1 = "SELECT flat_chemical, source_id \
                            FROM stitch.chemical_sources \
                            WHERE source_name = 'ATC' \
                            AND flat_chemical IN ('%s');" \
                            % "', '".join(list_of_stitch_ids)
                            
            mysql_query2 = "SELECT stereo_chemical, source_id \
                            FROM stitch.chemical_sources \
                            WHERE source_name = 'ATC' \
                            AND stereo_chemical IN ('%s');" \
                            % "', '".join(list_of_stitch_ids)
                            
            con = mdb.connect('localhost', 'root', 'kim630', 'stitch')
            with con:
                cur = con.cursor()
    
                cur.execute(mysql_query1)
                rows1 = cur.fetchall()
                
                cur.execute(mysql_query2)
                rows2 = cur.fetchall()
                
            self.drug_ATC = dict()
            for row in rows1:
                self.drug_ATC.setdefault(row[0], set()).add(row[1])
            for row in rows2:
                self.drug_ATC.setdefault(row[0], set()).add(row[1])
            
            with open('../data/pickled/drug_ATC.pickle', 'wb') as fh:
                pickle.dump(self.drug_ATC, fh)

            
        #######################################################################
        # Get the STITCH to Chembl mapping for the drugs in question
        if os.path.isfile('../data/pickled/list_of_chembl_ids.pickle'):
            with open('../data/pickled/list_of_chembl_ids.pickle', 'rb') as fh:
                list_of_chembl_ids = pickle.load(fh)
        else:
            list_of_chembl_ids = []
            for stitch_id in list_of_stitch_ids:
                # pubchem to chembl
                if stitch_id[:2] == 'CI': 
                    pubchem_id = str(int(stitch_id[4:]))
                    query_tuple = (pubchem_id, '22', '1')
                # drugbank to chembl
                elif stitch_id[:2] == 'DB':
                    pubchem_id = stitch_id
                    query_tuple = (pubchem_id, '2', '1')
                    
                try:
                    request_string = 'https://www.ebi.ac.uk/unichem/rest/src_compound_id_all/%s/%s/%s' % query_tuple
                    page = urllib2.urlopen(request_string).read()
                    unichem_output = json.loads(page)
                except urllib2.HTTPError as err:
                    print err
                    unichem_output = []
                    
                if unichem_output == []:
                    print stitch_id, pubchem_id, 'not found!'
                    list_of_chembl_ids.append('')
                else:
                    list_of_chembl_ids.append(str(unichem_output[0]['src_compound_id']))
                      
            with open('../data/pickled/list_of_chembl_ids.pickle', 'wb') as fh:
                pickle.dump(list_of_chembl_ids, fh)

        
        #######################################################################
        # Get the ATC codes using the Chembl database
        if os.path.isfile('../data/pickled/drug_ATC2.pickle'):
            with open('../data/pickled/drug_ATC2.pickle', 'rb') as fh:
                print "Number of ATC values before Chembl:", len(self.drug_ATC.values())
                self.drug_ATC = pickle.load(fh)
                print "Number of ATC values before Chembl:", len(self.drug_ATC.values())
                
        else:
            mysql_query = "SELECT md.chembl_id, mac.level5 \
                            FROM chembl_17.molecule_dictionary md \
                            INNER JOIN \
                            chembl_17.molecule_atc_classification mac \
                            ON md.molregno = mac.molregno \
                            WHERE md.chembl_id IN ('%s');" \
                            % "', '".join(list_of_chembl_ids)
                            
            con = mdb.connect('192.168.6.19', 'root', 'kim630', 'chembl_17')
            with con:
                cur = con.cursor()
                cur.execute(mysql_query)
                rows = cur.fetchall()
            
            for row in rows:
                stitch_id = list_of_stitch_ids[list_of_chembl_ids.index(row[0])]
                self.drug_ATC.setdefault(stitch_id, set()).add(row[1])
    
            with open('../data/pickled/drug_ATC2.pickle', 'wb') as fh:
                pickle.dump(self.drug_ATC, fh)
                
        all_drugs_info_df = pd.DataFrame({'stitch_id': list_of_stitch_ids,
                                  'pubchem_id': list_of_pubchem_ids,
                                  'chembl_id': list_of_chembl_ids,
                                  'SMILE': [ self.drug_SMILES.get(stitch_id, '') for stitch_id in list_of_stitch_ids ],
                                  'ATC': [ list(self.drug_ATC.get(stitch_id, '')) for stitch_id in list_of_stitch_ids ],
                                  'name': [ self.drug_names.get(stitch_id, '') for stitch_id in list_of_stitch_ids ]})
        with open('../all_drugs_info.tsv', 'w') as fh:
            all_drugs_info_df.to_csv(fh, sep='\t', 
                                     cols=['stitch_id', 'name', 'pubchem_id', 'chembl_id', 'SMILE', 'ATC'],
                                     index=False)
        
        
    
#    def load_saved_db_data(self):
#        """ Import internal data from separate pickled objects
#        """
#        with open('../data/pickled/' + 'drug_pairs' + '.pickle', 'rb') as ih:
#            self.drug_pairs = pickle.load(ih)
#        with open('../data/pickled/' + 'drugs_targets_pos' + '.pickle', 'rb') as ih:
#            self.drugs_targets_pos = pickle.load(ih)
#        with open('../data/pickled/' + 'drugs_targets_neg' + '.pickle', 'rb') as ih:
#            self.drugs_targets_neg = pickle.load(ih)
#        with open('../data/pickled/' + 'drug_data_pos' + '.pickle', 'rb') as ih:
#            self.drug_data_pos = pickle.load(ih)
#        with open('../data/pickled/' + 'target_data_pos' + '.pickle', 'rb') as ih:
#            self.target_data_pos = pickle.load(ih)
#        with open('../data/pickled/' + 'drug_data_neg' + '.pickle', 'rb') as ih:
#            self.drug_data_neg = pickle.load(ih)
#        with open('../data/pickled/' + 'target_data_neg' + '.pickle', 'rb') as ih:
#            self.target_data_neg = pickle.load(ih)
#        with open('../data/pickled/' + 'drug_SMILES' + '.pickle', 'rb') as ih:
#            self.drug_SMILES = pickle.load(ih)
#        
#
#    def save_db_data(self):
#        """ Export internal data as separate pickled objects
#        """
#        with open('../data/pickled/' + 'drug_pairs' + '.pickle', 'wb') as ih:
#            pickle.dump(self.drug_pairs, ih)
#        with open('../data/pickled/' + 'drugs_targets_pos' + '.pickle', 'wb') as ih:
#            pickle.dump(self.drugs_targets_pos, ih)
#        with open('../data/pickled/' + 'drugs_targets_neg' + '.pickle', 'wb') as ih:
#            pickle.dump(self.drugs_targets_neg, ih)
#        with open('../data/pickled/' + 'drug_data_pos' + '.pickle', 'wb') as ih:
#            pickle.dump(self.drug_data_pos, ih)
#        with open('../data/pickled/' + 'target_data_pos' + '.pickle', 'wb') as ih:
#            pickle.dump(self.target_data_pos, ih)
#        with open('../data/pickled/' + 'drug_data_neg' + '.pickle', 'wb') as ih:
#            pickle.dump(self.drug_data_neg, ih)
#        with open('../data/pickled/' + 'target_data_neg' + '.pickle', 'wb') as ih:
#            pickle.dump(self.target_data_neg, ih)
#        with open('../data/pickled/' + 'drug_SMILES' + '.pickle', 'wb') as ih:
#            pickle.dump(self.drug_SMILES, ih)


    
    def export_db_data_tsv(self):
        with open('../drug_data_pos.tsv', 'wb') as fh:
            self.drug_data_pos.to_csv(fh, sep='\t')
        with open('../drug_data_neg.tsv', 'wb') as fh:
            self.drug_data_neg.to_csv(fh, sep='\t')
            
        with open('../target_data_pos.tsv', 'wb') as fh:
            self.target_data_pos.to_csv(fh, sep='\t')
        with open('../target_data_neg.tsv', 'wb') as fh:
            self.target_data_neg.to_csv(fh, sep='\t')
            
        with open('../drug_data_merged_pos.tsv', 'wb') as fh:
            merged = p_db.drugs_targets_pos \
            .merge(p_db.drug_data_pos, how='left', left_on='drug_pair', right_index=True) \
            .merge(p_db.target_data_pos, how='left', left_on='target_pair', right_index=True)
            grouped = merged.groupby(['drug_pair'])
            maxed = grouped.max()
            maxed.to_csv(fh, sep='\t')
        with open('../drug_data_merged_neg.tsv', 'wb') as fh:
            merged = p_db.drugs_targets_neg \
            .merge(p_db.drug_data_neg, how='left', left_on='drug_pair', right_index=True) \
            .merge(p_db.target_data_neg, how='left', left_on='target_pair', right_index=True)
            grouped = merged.groupby(['drug_pair'])
            maxed = grouped.max()
            maxed.to_csv(fh, sep='\t')


    def _get_df(self, drug_or_target, pos_or_neg):
        """ Return a copy of the appropriate dataframe
        """
        if drug_or_target == 'drug' and pos_or_neg == 'pos':
            df_copy = self.drug_data_pos
        elif drug_or_target == 'drug' and pos_or_neg == 'neg':
            df_copy = self.drug_data_neg
        elif drug_or_target == 'target' and pos_or_neg == 'pos':
            df_copy = self.target_data_pos
        elif drug_or_target == 'target' and pos_or_neg == 'neg':
            df_copy = self.target_data_neg
        else:
            raise Exception("Wrong string entered for drug_or_target or pos_or_neg!")
            
        return df_copy

    
    def _set_df(self, drug_or_target, pos_or_neg, df_copy):
        """ Overwrite appropriate dataframe with a provided copy
        """
        if drug_or_target == 'drug' and pos_or_neg == 'pos':
            self.drug_data_pos = df_copy
        elif drug_or_target == 'drug' and pos_or_neg == 'neg':
            self.drug_data_neg = df_copy
        elif drug_or_target == 'target' and pos_or_neg == 'pos':
            self.target_data_pos = df_copy
        elif drug_or_target == 'target' and pos_or_neg == 'neg':
            self.target_data_neg = df_copy
        else:
            raise Exception("Wrong string entered for drug_or_target or pos_or_neg!")              
    
    
    def _read_data(self, filename, column_ids, feature_name):
        """
        Data is often supplied as separate files with two columns...
        The first column corresponds to the drug_pair or target_pair id, 
        and the second column corresponds to the feature values.
        This function reads such files into a list of ids and a list of values.
        """
        if feature_name in ['essentiality', 'shortest_path', 'atc_similarity', 'coexpression']:
            split_function = self._dummy_encode
        elif feature_name in ['pathway_shared_fraction']: 
            split_function = self._string_to_bool_encode_missing
        elif feature_name in ['functional_association', 'reg_network']:
            split_function = self._string_to_bool
        else:
            split_function = self._string_to_float
            
        pair_data = dict()
        with open(filename, 'r') as ih:
            for line in ih:
                row = line.split('\t')
                row1 = row[column_ids[0]].strip()
                row2 = row[column_ids[1]].strip()
                if row1 == '':
                    continue
                row2 = split_function(row2)
                # Take the maximum in  case of duplicates
                if pair_data.has_key(row1):
                    if any(row2 > pair_data[row1]):
                        pair_data[row1] = row2
                else:
                    pair_data[row1] = row2
                
        return pair_data.keys(), pair_data.values()


    def _string_to_float(self, value):
        try:
            result = float(value)
        except ValueError:
            return np.array([np.nan])
        return np.array([result,])


    def _string_to_bool(self, value):
        try:
            result = bool(float(value))
        except ValueError:
            raise Warning("String to bool converstion does not work with missing values!")
            return np.array([np.nan])
        # If some value
        if result:
            return np.array([1.])
        else:
            return np.array([0.])
            
    
    def _string_to_bool_encode_missing(self, value):
        try:
            result = bool(float(value))
        except ValueError:
            return np.array([0., 0.])
        # If some value
        if result:
            return np.array([1., 1.])
        else:
            return np.array([0., 1.])


    def _dummy_encode(self, value, num_cols=6):
        """ 
        One-encoding (using dummy variables) for gene essentiality. For each
        target, '-1' means non-essential, '0' means conflicting results, and '1'
        means essential. Assuming symetry for each target pair, there are 6 unique
        combinations.
        """
        feature = np.array([0.]*num_cols)
        
        try:
            ivalue = int(float(value))
        except ValueError:
            return feature

        if ivalue <= 1:
            feature[0] = 1
        elif ivalue == 2:
            feature[1] = 1
        elif ivalue == 3:
            feature[2] = 1           
        elif ivalue == 4:
            feature[3] = 1
        elif ivalue == 5:
            feature[4] = 1          
        elif ivalue >= 6:
            feature[5] = 1
        return feature

            
    def _gen_df_from_list(self, pair_index, pair_data, feature_name):
        """ 
        Add a list of values (with the list of corresponding indeces), to the
        dataframe that stores the data for all the features
        """        
        if len(np.shape(pair_data)) > 1:
            num_cols = np.shape(pair_data)[1]
        else:
            num_cols = 1
        
        # Get the column names for the new feature
        if num_cols == 1:
            feature_names = [ feature_name, ]
        elif num_cols > 1:
            feature_names = [ feature_name + '-' + str(i+1) for i in range(num_cols) ]
        
        # Create a dataframe for the new feature, and merge it with the old dataframe
        df_addition = pd.DataFrame(data=pair_data, index=pair_index, columns=feature_names)
        
        return df_addition
    
            
    def add_data(self, drug_or_target, pos_or_neg, filename, feature_name=None, 
                 na_value=None, random_pair=False, column_ids=[0,1]):
        """ 
        Add features to dataframes used for the internal storage of information
        """
        # Get a copy of the dataframe to which new a new feature will be added
        df_copy = self._get_df(drug_or_target, pos_or_neg)
        
        if not feature_name:
            # Read in a pandas dataframe that was saved in tsv format
            df_addition = pd.read_csv(filename, sep='\t', index_col=0)
            feature_name = df_addition.columns[0]
        else:
            # Get the index information (drug_pair or target_pair) and the feature data
            # for the new feature
            pair_index, pair_data = self._read_data(filename, column_ids, feature_name)
            assert len(pair_index) == len(pair_data)
            
            random.seed(131212)
            # Initially, each feature for the negative training set was calculated
            # using different random pairs of target proteins. If this is the case,
            # instead of using the supplied target_pair identifiers, generate a list
            # of identifiers randomly from the initialised dataframes. HACK!!!
            if random_pair:
                pair_index = list(df_copy.index)
                pair_data = pair_data * (len(pair_index)/len(pair_data) + 1)
                random.shuffle(pair_data) # Shuffle the list inplace
                pair_data = pair_data[:len(pair_index)]
            
            # If a feature is already present, don't add it again (this should not happen).
            if feature_name in [ column.split('-')[0] for column in df_copy.columns ]:
                raise Exception('Feature ' + feature_name + ' already exists!')
                
            # Add the features to the appropriate dataframe
            df_addition = self._gen_df_from_list(pair_index, pair_data, feature_name)


        new_columns = df_addition.columns
        df_copy = df_copy.merge(df_addition, how='left', left_index=True, right_index=True)
        
        for column in new_columns:
            print "%i added, %i Null, %i Zero, %i%%, to %s %s %s set" \
            % (len(df_copy[column]), sum(pd.isnull(df_copy[column])), sum(df_copy[column] == 0.),
            float(sum(pd.isnull(df_copy[column])) +  sum(df_copy[column] == 0.))/len(df_copy[column])*100,
            drug_or_target, column, pos_or_neg)
            
            if na_value == 'mean':
                df_copy[column].fillna(df_copy[column].mean(), inplace=True)
            elif na_value or na_value == 0.:
                df_copy[column].fillna(na_value, inplace=True)
                
        # Save the new dataframe
        self._set_df(drug_or_target, pos_or_neg, df_copy)


if __name__ == '__main__':
    
    p_db = project_db()

    if os.path.isfile('../data/pickled/ALLDB.pickle'):
        with open('../data/pickled/ALLDB.pickle', 'r') as fh:
            p_db = pickle.load(fh)
            p_db.export_db_data_tsv()
            
    else:        
        ###########################################################################
        # Data for drug pairs
        
        p_db.add_data('drug', 'pos', '../data/drugs/chemical-similarity/chemical_similarity_data_pos.tsv')
        p_db.add_data('drug', 'neg', '../data/drugs/chemical-similarity/chemical_similarity_data_neg.tsv')
        
        p_db.add_data('drug', 'pos', '../data/drugs/clinical-similarity/atc_similarity_pos.tsv', 'atc_similarity', na_value=0) # dummy encoding
        p_db.add_data('drug', 'neg', '../data/drugs/clinical-similarity/atc_similarity_neg.tsv', 'atc_similarity', na_value=0) # dummy encoding
    
        
        ###########################################################################
        # Date for target pairs
        
        p_db.add_data('target', 'pos', '../data/targets/Coexpression/coexpression_pos.tsv', 'coexpression', na_value=0.) # dummy encoding
        p_db.add_data('target', 'neg', '../data/targets/Coexpression/coexpression_neg.tsv', 'coexpression', na_value=0.) # dummy encoding
        
        p_db.add_data('target', 'pos', '../data/targets/FxnAssociation/protein.links.v9.05_Human_ENSP_ENSG.txt', 'functional_association', na_value=0., column_ids=[0,2]) # bool encoding
        p_db.add_data('target', 'neg', '../data/targets/FxnAssociation/protein.links.v9.05_Human_ENSP_ENSG.txt', 'functional_association', na_value=0., column_ids=[0,2]) # bool encoding
        
        p_db.add_data('target', 'pos', '../data/targets/Essentiality/essentiality_pos.tsv', 'essentiality', na_value=0.) # dummy encoding
        p_db.add_data('target', 'neg', '../data/targets/Essentiality/essentiality_neg.tsv', 'essentiality', na_value=0.) # dummy encoding
        
        p_db.add_data('target', 'pos', '../data/targets/ShortestPath/CRG_ShortestPathLength-pos.txt', 'shortest_path', na_value=0.) # dummy encoding
        p_db.add_data('target', 'neg', '../data/targets/ShortestPath/CRG_ShortestPathLength-neg.txt', 'shortest_path', na_value=0.) # dummy encoding
    
        
        ###########################################################################
        # INCORRECT TRAINING POINTS FROM THIS POINT ON!!!
        
        p_db.add_data('target', 'pos', '../data/targets/Pathway_SharedFraction/Pathway_SharedFraction-pos.txt', 'pathway_shared_fraction', na_value=0.) # bool encoding
        p_db.add_data('target', 'neg', '../data/targets/Pathway_SharedFraction/Pathway_SharedFraction-neg.txt', 'pathway_shared_fraction', na_value=0., random_pair=True)  # bool encoding
        
        p_db.add_data('target', 'pos', '../data/targets/GOAll_SharedFraction/GOAll_SharedFraction-pos.txt', 'go_all_shared_fraction')
        p_db.add_data('target', 'neg', '../data/targets/GOAll_SharedFraction/GOAll_SharedFraction-neg.txt', 'go_all_shared_fraction', random_pair=True)
        
    #    p_db.add_data('target', 'pos', data_dir + '/targets/GOBP_Pattern/GO_BP_Pattern-pos.txt', 'go_bp_pattern')
    #    p_db.add_data('target', 'neg', data_dir + '/targets/GOBP_Pattern/GO_BP_Pattern-neg.txt', 'go_bp_pattern', random_pair=True)
        
        p_db.add_data('target', 'pos', '../data/targets/GOBP_SharedFraction/GOBP_SharedFraction-pos.txt', 'go_bp_shared_fraction')
        p_db.add_data('target', 'neg', '../data/targets/GOBP_SharedFraction/GOBP_SharedFraction-neg.txt', 'go_bp_shared_fraction', random_pair=True)
        
    #    p_db.add_data('target', 'pos', data_dir + '/targets/GOCC_Pattern/GO_CC_Pattern-pos.txt', 'go_cc_pattern')
    #    p_db.add_data('target', 'neg', data_dir + '/targets/GOCC_Pattern/GO_CC_Pattern-neg.txt', 'go_cc_pattern', random_pair=True)
        
        p_db.add_data('target', 'pos', '../data/targets/GOCC_SharedFraction/GOCC_SharedFraction-pos.txt', 'go_cc_shared_fraction')
        p_db.add_data('target', 'neg', '../data/targets/GOCC_SharedFraction/GOCC_SharedFraction-neg.txt', 'go_cc_shared_fraction', random_pair=True)
        
    #    p_db.add_data('target', 'pos', data_dir + '/targets/GOMF_Pattern/GO_MF_Pattern-pos.txt', 'go_mf_pattern')
    #    p_db.add_data('target', 'neg', data_dir + '/targets/GOMF_Pattern/GO_MF_Pattern-neg.txt', 'go_mf_pattern', random_pair=True)
        
        p_db.add_data('target', 'pos', '../data/targets/GOMF_SharedFraction/GOMF_SharedFraction-pos.txt', 'go_mf_shared_fraction')
        p_db.add_data('target', 'neg', '../data/targets/GOMF_SharedFraction/GOMF_SharedFraction-neg.txt', 'go_mf_shared_fraction', random_pair=True)
        
    #    p_db.add_data('target', 'pos', data_dir + '/targets/Pathway_Pattern/Pathway_Pattern-pos.txt', 'pathway_pattern')
    #    p_db.add_data('target', 'neg', data_dir + '/targets/Pathway_Pattern/Pathway_Pattern-neg.txt', 'pathway_pattern', random_pair=True)  
        
        p_db.add_data('target', 'pos', '../data/targets/RegNetwork/RegNetwork-pos.txt', 'reg_network', na_value=0.) # bool encoding
        p_db.add_data('target', 'neg', '../data/targets/RegNetwork/RegNetwork-neg.txt', 'reg_network', na_value=0., random_pair=True) # bool encoding
        
        ###########################################################################
        # TOO FEW DATA POINTS
        
    #    p_db.add_data('target', 'pos', data_dir + '/targets/GeneticInteraction/GeneticInteraction-pos.txt', 'genetic_interaction')
    #    p_db.add_data('target', 'neg', data_dir + '/targets/GeneticInteraction/GeneticInteraction-neg.txt', 'genetic_interaction', random_pair=True)
        
        p_db.export_db_data_tsv()
        with open('../data/pickled/ALLDB.pickle', 'w') as fh:
            pickle.dump(p_db, fh)