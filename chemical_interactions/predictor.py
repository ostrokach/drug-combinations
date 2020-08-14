# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 12:27:21 2014

@author: alexey
"""

import os
import logging
import numpy as np
import pandas as pd
from collections import deque, defaultdict, OrderedDict
from sklearn.ensemble import GradientBoostingClassifier
import imp



#%% Parameters
#2014-04-17 best parameters
predictor_parameters_all = {
    'predictor_1': {
        'clf_name': 'Predictor 1',
        'clf_type': 1,
        'clf_idx': 0,
        'plot_col_idx': 0,
        'input_file': 'predictor_1.tsv',
        'clf_options': {
            'loss': 'deviance',
            'n_estimators': 900,
            'learning_rate': 0.01,
            'max_depth': 2,
            'subsample': 0.3}},

    'predictor_1hs': {
        'clf_name': 'Predictor 1 High Scored',
        'clf_type': 1,
        'clf_idx': 0,
        'plot_col_idx': 1,
        'input_file': 'predictor_1_high_scored.tsv',
        'clf_options': {
            'loss': 'deviance',
            'n_estimators': 700,
            'learning_rate': 0.05,
            'max_depth': 2,
            'subsample': 0.5}},

    'predictor_2': {
        'clf_name': 'Predictor 2',
        'clf_type': 2,
        'clf_idx': 5,
        'plot_col_idx': 2,
        'input_file': 'predictor_2.tsv',
        'clf_options': {
            'loss': 'deviance',
            'n_estimators': 900,
            'learning_rate': 0.01,
            'max_depth': 4,
            'subsample': 0.7}},

    'predictor_2hs': {
        'clf_name': 'Predictor 2 High Scored',
        'clf_type': 2,
        'clf_idx': 5,
        'plot_col_idx': 3,
        'input_file': 'predictor_2_high_scored.tsv',
        'clf_options': {
            'loss': 'deviance',
            'n_estimators': 700,
            'learning_rate': 0.1,
            'max_depth': 4,
            'subsample': 0.7}}
}


columnname_mapping = OrderedDict([
    ('Type', 'Type'),
    ('DrugPair', 'DrugPair'),
    ('TargetPair', 'TargetPair'),
    ('labels', 'Labels'),
    ('y_true_all', 'y_true_all'),
    ('y_pred_all', 'y_pred_all'),
    ('probas_p1', 'probas_p1'),
    ('probas_p1hs', 'probas_p1hs'),
    ('probas_p2', 'probas_p2'),
    ('probas_p2hs', 'probas_p2hs'),
    ('morganfingerprintr2_tanimoto', 'chemical_similarity'),
    ('atc_similarity', 'atc_similarity'),
    ('side_effect_similarity', 'side_effect_similarity'),
    ('coexpression', 'gene_coexpression'),
    ('essentiality', 'gene_coessentiality'),
    ('shortest_path_length', 'biogrid_shortest_path_length'),
    ('eb_max', 'biogrid_eb_max'),
    ('shortest_path_length.1', 'getint_shortest_path_length'),
    ('eb_max.1', 'getint_eb_max'),
    ('shortest_path_length.2', 'string_shortest_path_length'),
    ('eb_max.2', 'string_eb_max'),
    ('gene_essentiality_1', 'gene_essentiality_1'),
    ('gene_essentiality_2', 'gene_essentiality_2'),
    ('go_all_sem_sim', 'go_all_sem_sim'),
    ('go_bp_sem_sim', 'go_bp_sem_sim'),
    ('go_cc_sem_sim', 'go_cc_sem_sim'),
    ('go_mf_sem_sim', 'go_mf_sem_sim'),
    ('phylogenic_similarity', 'phylogenic_similarity'),
])



#%% Functions
def get_confirmed_testfile(input_testfile, confirmed_drugpairs, input_folder, suffix=''):
    input_file_data = pd.read_csv(os.path.join(input_folder, input_testfile), sep='\t')
    input_file_data.columns = [columnname_mapping.get(column, column) for column in input_file_data.columns]

    curated_drug_pairs = pd.read_csv(input_folder + confirmed_drugpairs, sep='\t')
    curated_drug_pairs['cid_1'] = [int(drug_pair.split('_')[0][3:]) for drug_pair in curated_drug_pairs.DrugPair.values]
    curated_drug_pairs['cid_2'] = [int(drug_pair.split('_')[1][3:]) for drug_pair in curated_drug_pairs.DrugPair.values]
    assert sum([x[0] > x[1] for x in curated_drug_pairs[['cid_1', 'cid_2']].values]) == 0
    curated_drug_pairs = curated_drug_pairs[['cid_1','cid_2']]

    input_file_data = input_file_data.merge(curated_drug_pairs, on=['cid_1', 'cid_2'])
    confirmed_testfile = input_testfile + '.confirmed' + suffix
    input_file_data.to_csv(input_folder + confirmed_testfile, sep='\t', index=False)
    return confirmed_testfile



#%% Classes
class Predictor(object):

    def __init__(self, input_file, input_folder, output_folder, random_state=None, logger=None):
        """ Initialize the predictor
        """
        if random_state is None:
            self.random_state = np.random.RandomState(42)
        else:
            self.random_state = random_state

        if logger is None:
            imp.reload(logging)
            self.logger = logging.getLogger("chemical_interactions")
            self.logger.handlers = [logging.StreamHandler()]
            self.logger.setLevel(logging.DEBUG)
            self.logger.debug('Log successfully initialized')
        else:
            self.logger = logger

        self.input_file = input_file
        self.input_folder = input_folder
        self.output_folder = output_folder

        input_file = os.path.join(self.input_folder, self.input_file)
        self.logger.info('Reading data from file: {}'.format(input_file))
        predictor_df = pd.read_csv(input_file, sep='\t')

        # Fill missing values for GO annotations with 0
        ###predictor_df.fillna(0.0, inplace=True)
        predictor_df.loc[:,'go_bp_sem_sim'].fillna(0.0, inplace=True)
        predictor_df.loc[:,'go_cc_sem_sim'].fillna(0.0, inplace=True)
        predictor_df.loc[:,'go_mf_sem_sim'].fillna(0.0, inplace=True)
        predictor_df.loc[:,'go_all_sem_sim'].fillna(0.0, inplace=True)
        # string_topo_eb
        # Not really sure if ok to add gene essentiality, but ok...
        predictor_df.loc[:,'gene_essentiality_1'].fillna(0, inplace=True)
        predictor_df.loc[:,'gene_essentiality_2'].fillna(0, inplace=True)

        # Shuffle the dataframe and then sort it by `Type` so that x-fold cross-validation is done
        # using a random subset with an equal number of positive and negative training examples.
        if 'Type' in predictor_df.columns:
            predictor_df = predictor_df.sort(['Type'])

        # Limit the number of rows that can be considered for the training set
        #if len(predictor_df) > int(1e5):
        #   new_index = self.random_state.permutation(predictor_df.index)
        #   predictor_df = predictor_df.iloc[new_index[:int(1e5)]]

        # Add `TargetPair` and `DrugPair` columns which are absent from data files extracted from
        # the database
        self.logger.debug('Adding TargetPair column...')
        predictor_df['TargetPair'] = [
            '_'.join(['ENSP' + str(x) for x in xx])
            for xx
            in predictor_df[['ensp_1', 'ensp_2']].values]
        del predictor_df['ensp_1']
        del predictor_df['ensp_2']

        if 'cid_1' in predictor_df.columns and 'cid_2' in predictor_df.columns:
            self.logger.debug('Adding DrugPair column...')
            predictor_df['DrugPair'] = [
                '_'.join(['CID' + str(x) for x in xx])
                for xx
                in predictor_df[['cid_1', 'cid_2']].values]
            del predictor_df['cid_1']
            del predictor_df['cid_2']

        # Format the `gene essentiality feature`
        self.logger.debug('Reformatting features...')
        predictor_df['essentiality'] = [
            '_'.join([str(int(x)) for x in xx])
            for xx
            in predictor_df[['gene_essentiality_1', 'gene_essentiality_2']].values]
        predictor_df['essentiality'] = predictor_df['essentiality'].apply(self.categorize_essentiality)
        del predictor_df['gene_essentiality_1']
        del predictor_df['gene_essentiality_2']

        self.predictor_df = predictor_df


    @staticmethod
    def categorize_essentiality(feature):
        """ Assign an integer value to each possible gene essentiality combinations
        """
        feature = [int(i) for i in feature.split('_')]
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


    def cross_validate_predictor(self, clf_options, n_folds=False, columns_to_keep=[], columns_to_drop=[], additional_columns_to_drop=[]):
        """
        Cross-validate the predictor, and return predicted and actual scores for all points
        in the data set. In addition, return a predictor trained using all availible data.
        """
        df = self.predictor_df

        # Obtain indices for leave-one-out cross-validation of the predictor
        # This may take a while, so store the indices for future retreival
        xval_1fold_indices_filename = self.output_folder + self.input_file + '.1fold_xval.npy'
        if os.path.isfile(xval_1fold_indices_filename):
            self.pos_neg_indices = np.load(xval_1fold_indices_filename)
        else:
            self.pos_neg_indices = self.leave_drug_pair_out()
            np.save(xval_1fold_indices_filename, self.pos_neg_indices)

        # If doing x-fold cross-validation, obtain index arrays for a given `x-fold`
        if n_folds and n_folds > 1:
            xval_kfold_indices_filename = self.output_folder + self.input_file + '.%ifold_xval.npy' % n_folds
            if os.path.isfile(xval_kfold_indices_filename):
                self.pos_neg_indices_kfold = np.load(xval_kfold_indices_filename)
            else:
                self.pos_neg_indices_kfold = self.join_drug_pair_indices(self.pos_neg_indices, n_folds)
                np.save(xval_kfold_indices_filename, self.pos_neg_indices_kfold)
            kf_xv = self.pos_neg_indices_kfold
        else:
            kf_xv = self.pos_neg_indices

        # Get names of columns that have features to be used by the predictor
        if not columns_to_keep and columns_to_drop:
            columns_to_keep = [col for col in df.columns if col not in columns_to_drop]
        elif not columns_to_keep and not columns_to_drop:
            columns_to_drop = [
                'Type',
                'TargetPair',
                'DrugPair'] + additional_columns_to_drop
            columns_to_keep = [col for col in df.columns if col not in columns_to_drop]
        self.logger.debug('`columns_to_keep`:\n' + ' '.join(sorted(columns_to_keep)))
        self.logger.debug('`predictor_df.columns`:\n' + ' '.join(sorted(df.columns)))

        # Get the data and labels as numpy arrays
        data, labels = self.get_x_y(df, columns_to_keep)

        # Pre-allocate bjects for holding cross-validation results
        y_pred_all = []
        y_true_all = []
        dfs_left_out = []
        preds_per_drug_pair = defaultdict(list)

        # Create a classifier
        clf_options.update({'random_state': self.random_state})
        clf = GradientBoostingClassifier(**clf_options)

        # Perform cross-validation
        for xv_idx, (test, train) in enumerate(kf_xv):
            # Skip cases where training set does not have both positive and negative
            # examples.
            if not (any(labels[train]==True) and any(labels[train]==False)):
                self.logger.debug('skipping no. %i' % xv_idx)
                continue

            clf.fit(data[train,:], labels[train])

            probas_ = clf.predict_proba(data[test,:])
            y_pred_all.extend(probas_[:,1])
            y_true_all.extend(labels[test])
            df_left_out = df.iloc[np.where(test)]
            df_left_out['y_pred_all'] = probas_[:,1]
            df_left_out['y_true_all'] = labels[test]
            dfs_left_out.append(df_left_out)

            # Combine scores for all target pairs for a given drug pair
            # (for `y_true_all_perdrugpair`, `y_pred_all_perdrugpair`)
            df_row_nums = []
            for df_row_num, test_idx in enumerate(test):
                if test_idx:
                    df_row_nums.append(df_row_num)
            if 'DrugPair' in df.columns:
                for df_row_num, proba in zip(df_row_nums, probas_[:,1]):
                    drug_pair = df.iloc[df_row_num]['DrugPair']
                    preds_per_drug_pair[drug_pair].append(proba)

            self.logger.info('Cross-validating round %i out of %i' % (xv_idx+1, len(kf_xv)))

        # Combine predictions for all data points
        y_pred_all = np.array(y_pred_all, dtype=float)
        y_true_all = np.array(y_true_all, dtype=bool)
        dfs_left_out = pd.concat(dfs_left_out, ignore_index=True)

        # Compute the maximum over all target pairs for each drug pair
        if 'DrugPair' in df.columns:
            temp = [ (list(df[df['DrugPair']==x[0]]['Type'])[0]=='Pos', max(x[1]),) for x in list(preds_per_drug_pair.items()) ]
            temp = np.array(temp)
            y_true_all_perdrugpair = temp[:,0]
            y_pred_all_perdrugpair = temp[:,1]
        else:
            y_true_all_perdrugpair = None
            y_pred_all_perdrugpair = None

        clf.fit(data, labels)
        return clf, y_true_all, y_pred_all, y_true_all_perdrugpair, y_pred_all_perdrugpair, dfs_left_out


    def get_x_y(self, df, columns_to_keep):
        x = np.float64(df[columns_to_keep].values)
        y = df['Type'] == 'Pos'
        return x, y


    def leave_drug_pair_out(self):
        df = self.predictor_df
        global temp
        pos_neg_indices = []
        temp = pos_neg_indices
        # Get positive and negative indices for each drug pair
        if 'DrugPair' in self.predictor_df.columns:
            self.logger.debug('Have a DrugPair column. Doing leave one drug pair out xval...')
            visited_drug_pairs = set()
            for drug_1, drug_2 in (drug_pair.split('_') for drug_pair in df['DrugPair'].values):
                if (drug_1, drug_2) in visited_drug_pairs:
                    continue
                visited_drug_pairs.add((drug_1, drug_2))

                index_pos = [ (drug_1 in drug_pair.split('_')) and (drug_2 in drug_pair.split('_'))
                    for drug_pair in df['DrugPair'].values ]
                drugs_to_leave_out = set()
                targets_to_leave_out = set()
                deque((drugs_to_leave_out.update(drug_pair.split('_')) for drug_pair in df['DrugPair'].iloc[index_pos]), maxlen=0)
                deque((targets_to_leave_out.update(target_pair.split('_')) for target_pair in df['TargetPair'].iloc[index_pos]), maxlen=0)
                index_neg = [
                    not ((drug_pair.split('_')[0] in drugs_to_leave_out)
                    or (drug_pair.split('_')[1] in drugs_to_leave_out)
                    or (target_pair.split('_')[0] in targets_to_leave_out)
                    or (target_pair.split('_')[1] in targets_to_leave_out))
                    for (drug_pair, target_pair) in zip(df['DrugPair'], df['TargetPair'])]
                pos_neg_indices.append([np.array(index_pos, dtype=bool), np.array(index_neg, dtype=bool)])
        else:
            self.logger.debug('Don\'t have a DrugPair column! Doing leave one target pair out xval...')
            for target_1, target_2 in (target_pair.split('_') for target_pair in df['TargetPair'].values):
                index_pos = [ (target_1 in target_pair.split('_')) and (target_2 in target_pair.split('_'))
                    for target_pair in df['TargetPair'].values ]
                index_neg = [ (target_1 not in target_pair.split('_')) and (target_2 not in target_pair.split('_'))
                    for target_pair in df['TargetPair'].values ]
                pos_neg_indices.append([np.array(index_pos, dtype=bool), np.array(index_neg, dtype=bool)])
        return pos_neg_indices


    def join_drug_pair_indices(self, pos_neg_indices, n_folds):
        items_per_fold = int(np.ceil(float(len(pos_neg_indices))/n_folds))
        n_folds = int(np.ceil(float(len(pos_neg_indices))/items_per_fold))
        pos_neg_indices_xval = [[]] * n_folds
        for i in range(items_per_fold):
            for j in range(n_folds):
                if j+i*n_folds >= len(pos_neg_indices):
                    continue
                elif i == 0:
                    index_pos, index_neg = pos_neg_indices[j]
                    pos_neg_indices_xval[j] = [index_pos, index_neg]
                else:
                    index_pos_1, index_neg_1 = pos_neg_indices_xval[j]
                    index_pos_2, index_neg_2 = pos_neg_indices[j+i*n_folds]
                    index_pos = index_pos_1 | index_pos_2
                    index_neg = index_neg_1 & index_neg_2
                    pos_neg_indices_xval[j] = [index_pos, index_neg]
        return pos_neg_indices_xval


    def get_data_and_labels(self, additional_columns_to_drop=[]):
        columns_to_drop = ['Type', 'TargetPair', 'DrugPair'] + additional_columns_to_drop
        columns_to_keep = [col for col in self.predictor_df.columns if col not in columns_to_drop]
        data = np.float64(self.predictor_df[columns_to_keep].values)
        if 'Type' in self.predictor_df.columns:
            labels = self.predictor_df['Type'] == 'Pos'
        else:
            labels = np.zeros( (data.shape[0],1), dtype=bool )
        return data, labels
