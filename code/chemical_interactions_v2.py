# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 12:27:21 2014

@author: alexey
"""

import os
import numpy as np
import scipy as sp
import pandas as pd
import math

from collections import defaultdict

pd.options.mode.chained_assignment = None


#%% Initiate the logger
import logging
logger = logging.getLogger("chemical_interactions")
logger.handlers = []
logger.setLevel(logging.DEBUG)
logger.addHandler(logging.StreamHandler())
logger.debug('log successfully initialized')
random_state = np.random.RandomState(42)



#%% Parameters for external use. Do not remove!
n_folds = 40
path_to_data = '/home/kimlab1/strokach/databases/chemical_interactions/'
output_folder = 'version_2.1/final_predictors_2/'
input_files = [
    'version_2.1/predictors_2/predictor_1.tsv',
    'version_2.1/predictors_2/predictor_1_high_scored.tsv',
    'version_2.1/predictors_2/predictor_1_independent_validation.tsv',
    'version_2.1/predictors_2/predictor_1_all_unused.tsv',
    'version_2.1/predictors_2/predictor_1_all_unused_pairs.tsv',
    'version_2.1/predictors_2/predictor_1_all_unused_pairs_3.tsv',

    'version_2.1/predictors_2/predictor_2.tsv',
    'version_2.1/predictors_2/predictor_2_high_scored.tsv',
    'version_2.1/predictors_2/predictor_2_independent_validation.tsv',
    'version_2.1/predictors_2/predictor_2_all_unused.tsv',
    'version_2.1/predictors_2/predictor_2_all_unused_pairs.tsv',
    'version_2.1/predictors_2/predictor_2_all_unused_pairs_3.tsv',   

    'version_2.1/predictors_2/predictor_0.tsv',
    '', #'version_2.1/predictors_2/predictor_0_high_scored.tsv',
    'version_2.1/predictors_2/predictor_0_independent_validation.tsv',
    '', #'version_2.1/predictors_2/predictor_0_all_unused.tsv',
    'version_2.1/predictors_2/predictor_0_all_unused_pairs.tsv',
    'version_2.1/predictors_2/predictor_0_all_unused_pairs_3.tsv'
    ]



#%% Column name variables
header_columns = [
    'Type',
    'TargetPair',
    'DrugPair']

chemical_similarity_columns = [
    'RDKFingerprint_Tanimoto', #0
    'RDKFingerprint_Dice',
    'RDKFingerprint_Cosine',
    'RDKFingerprint_Russel',
    'RDKFingerprint_Kulczynski',
    'RDKFingerprint_McConnaughey', #5
    'FingerprintMol_Tanimoto',
    'FingerprintMol_Dice',
    'FingerprintMol_Cosine',
    'FingerprintMol_Russel',
    'FingerprintMol_Kulczynski', #10
    'FingerprintMol_McConnaughey',
    'MACCSkeys_Tanimoto',
    'MACCSkeys_Dice',
    'MACCSkeys_Cosine',
    'MACCSkeys_Russel', #15
    'MACCSkeys_Kulczynski',
    'MACCSkeys_McConnaughey',
    'AtomPairFingerprint_Tanimoto',
    'AtomPairFingerprint_Dice',
    'AtomPairFingerprint_Cosine', #20
    'AtomPairFingerprint_Russel',
    'AtomPairFingerprint_Kulczynski',
    'AtomPairFingerprint_McConnaughey',
    'TopologicalTorsionFingerprint_Tanimoto',
    'TopologicalTorsionFingerprint_Dice', #25
    'MorganFingerprintR2_Tanimoto',
    'MorganFingerprintR2_Dice',
    'MorganFingerprintR2withFeatures_Tanimoto',
    'MorganFingerprintR2withFeatures_Dice'] #29

chemical_other_columns = [
    'atc_similarity',
    'side_effect_similarity']

target_pairwise_columns = [
    'ShortestPathLength',
    'Essentiality',
    'Coexpression',
    'ShortestPathLength_GetInt',
    'All',
    'BP',
    'MF',
    'CC',
    'Phylo',
    'ShortestPathLength_PPI',
    'EB_Maximum',
    'EB_Maximum_GetInt_EB',
    'EB_Maximum_STRINGTopo_EB']

target_pairwise_columns_supl = [
    'Slim',
    'EB_Minimum',
    'EB_Mean',
    'EB_Fraction',
    'EB_Minimum_GetInt_EB',
    'EB_Mean_GetInt_EB',
    'EB_Fraction_GetInt_EB',
    'EB_Minimum_STRINGTopo_EB',
    'EB_Mean_STRINGTopo_EB',
    'EB_Fraction_STRINGTopo_EB']

target_nonpairwise_columns = [
    'Degree',
    'ClusteringCoeff',
    'Betweenness',
    'Closeness',
    'Degree_GetInt',
    'ClusteringCoeff_GetInt',
    'Betweenness_GetInt',
    'Closeness_GetInt',
    'Degree_PPI',
    'ClusteringCoeff_PPI',
    'Betweenness_PPI',
    'Closeness_PPI']



###############################################################################
#%% Import data from all files into a single pandas dataframe and run ML

class Predictor(object):

    def __init__(self, input_file, path_to_data, clf=None):

        logger.debug(input_file)
        self.path_to_data = path_to_data
        self.input_file = input_file
        self.clf = clf
        self.predictor_df = pd.read_csv(path_to_data + input_file, sep='\t')

        # Take only the first 100 thousand columns if it's a big file
#        if len(self.predictor_df) > 1e6:
#            new_index = random_state.permutation(self.predictor_df.index)
#            self.predictor_df = self.predictor_df.iloc[new_index[:1e5]]

        self.predictor_df.fillna(0.0, inplace=True)

        if 'Type' in self.predictor_df.columns:
            self.predictor_df = self.predictor_df.sort(['Type'])

        if (('ensp_1' in self.predictor_df.columns) and
            ('ensp_2' in self.predictor_df.columns)):
                self.predictor_df['TargetPair'] = self.predictor_df[['ensp_1', 'ensp_2']].apply(lambda x: 'ENSP{}_ENSP{}'.format(*x), axis=1)
                del self.predictor_df['ensp_1']
                del self.predictor_df['ensp_2']

        if (('cid_1' in self.predictor_df.columns) and
            ('cid_2' in self.predictor_df.columns)):
                logger.debug('Adding DrugPair column...')
                self.predictor_df['DrugPair'] = self.predictor_df[['cid_1', 'cid_2']].apply(lambda x: 'CID{}_CID{}'.format(*x), axis=1)
                del self.predictor_df['cid_1']
                del self.predictor_df['cid_2']

        if (('gene_essentiality_1' in self.predictor_df.columns) and
            ('gene_essentiality_2' in self.predictor_df.columns)):
                self.predictor_df['essentiality'] = self.predictor_df[['gene_essentiality_1', 'gene_essentiality_2']].apply(lambda x: '{}_{}'.format(*x), axis=1)
                self.predictor_df['essentiality'] = self.predictor_df['essentiality'].apply(self._categorize_essentiality)
                del self.predictor_df['gene_essentiality_1']
                del self.predictor_df['gene_essentiality_2']


    def _categorize_essentiality(self, feature):
        """ Assign an integer on the scale of 1 to 10 to describe all possible
        combinations of gene essentiality scores. """
        try:
            feature = [int(i) for i in feature.split('_')]
        except AttributeError:
            logger.debug(feature)
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


    def get_data_and_labels(self):
        columns_to_drop = ['Type', 'TargetPair', 'DrugPair']
        columns_to_keep = [col for col in self.predictor_df.columns if col not in columns_to_drop]
        data = np.float64(self.predictor_df[columns_to_keep].values)
        if 'Type' in self.predictor_df.columns:
            labels = self.predictor_df['Type'] == 'Pos'
        else:
            labels = np.zeros( (data.shape[0],1), dtype=bool )
        return data, labels


    def cross_validate_predictor(self, n_folds=None, columns_to_keep=[], columns_to_drop=[], additional_columns_to_drop=[]):
        """
        """
        # Create positive / negative cross-validation splits
        xval_1fold_indices_filename = self.path_to_data + self.input_file + '.1fold_xval.npy'
        if os.path.isfile(xval_1fold_indices_filename):
            self.pos_neg_indices = np.load(xval_1fold_indices_filename)
        else:
            self.pos_neg_indices = self._leave_drug_pair_out()
            np.save(xval_1fold_indices_filename, self.pos_neg_indices)

        if n_folds and n_folds > 1:
            xval_kfold_indices_filename = self.path_to_data + self.input_file + '.%ifold_xval.npy' % n_folds
            if os.path.isfile(xval_kfold_indices_filename):
                self.pos_neg_indices_kfold = np.load(xval_kfold_indices_filename)
            else:
                self.pos_neg_indices_kfold = self.join_drug_pair_indices(self.pos_neg_indices, n_folds)
                np.save(xval_kfold_indices_filename, self.pos_neg_indices_kfold)
            kf_xv = self.pos_neg_indices_kfold
        else:
            kf_xv = self.pos_neg_indices

        # Get a list of columns containing features to be used by the predicor
        if not columns_to_keep and columns_to_drop:
            columns_to_keep = [col for col in self.predictor_df.columns if col not in columns_to_drop]
        elif not columns_to_keep and not columns_to_drop:
            columns_to_drop = [
                u'Type',
                u'TargetPair',
                u'DrugPair'] + additional_columns_to_drop
            columns_to_keep = [col for col in self.predictor_df.columns if col not in columns_to_drop]
        logger.debug(self.predictor_df.columns)
        logger.debug(columns_to_keep)

        # Get data and labels numpy arrays
        data, labels = self.get_x_y(self.predictor_df, columns_to_keep)

        # Keep the predicted and actual values for all points in the validation set
        y_pred_all = []
        y_true_all = []
        dfs_left_out = []

        # Keep the number of samples in the training and validation sets in order
        # to track positive / negative set inbalance
        number_of_samples_in_positive_and_negative_training_sets = []
        number_of_samples_in_positive_and_negative_validation_sets = []
        preds_per_drug_pair = defaultdict(list)

        # Loop through teach train/test split in the cross-validation matrix
        for xv_idx, (test, train) in enumerate(kf_xv):
            if (not any(labels[train]==True) or
                not any(labels[train]==False)):
                    logger.debug('skipping no. %i' % xv_idx)
                    continue

            do_weighted = False
            if do_weighted:
                sample_weight = np.ones(len(labels[train]))
                sample_weight[np.array(labels[train]==1, dtype=bool)] *= float(sum(labels[train]==0))/float(sum(labels[train]==1))
                logger.debug(sum(sample_weight[np.array(labels[train]==1, dtype=bool)]))
                logger.debug(sum(sample_weight[np.array(labels[train]==0, dtype=bool)]))
                self.clf.fit(data[train,:], labels[train], sample_weight=sample_weight)
            else:
                self.clf.fit(data[train,:], labels[train])

            probas_ = self.clf.predict_proba(data[test,:])
            y_pred_all.extend(probas_[:,1])
            y_true_all.extend(labels[test])
            df_left_out = self.predictor_df.iloc[np.where(test)]
            df_left_out['y_pred_all'] = probas_[:,1]
            df_left_out['y_true_all'] = labels[test]
            dfs_left_out.append(df_left_out)

            df_row_nums = []
            for df_row_num, test_idx in enumerate(test):
                if test_idx:
                    df_row_nums.append(df_row_num)
            if 'DrugPair' in self.predictor_df.columns:
                for df_row_num, proba in zip(df_row_nums, probas_[:,1]):
                    drug_pair = self.predictor_df.iloc[df_row_num]['DrugPair']
                    preds_per_drug_pair[drug_pair].append(proba)

            logger.info('Cross-validating round %i out of %i' % (xv_idx+1, len(kf_xv)))
            number_of_samples_in_positive_and_negative_training_sets.append([sum(labels[train]==0), sum(labels[train]==1)])
            number_of_samples_in_positive_and_negative_validation_sets.append([sum(labels[test]==0), sum(labels[test]==1)])


        # Print some statistics
        logger.info('Mean number of samples in the training set (pos/neg): ({}/{})'.format(
            sum(zip(*number_of_samples_in_positive_and_negative_training_sets)[0]) /
                float(len(number_of_samples_in_positive_and_negative_training_sets)),
            sum(zip(*number_of_samples_in_positive_and_negative_training_sets)[1]) /
                float(len(number_of_samples_in_positive_and_negative_training_sets))))

        logger.info('Number of samples in test set (pos/neg): ({}/{})'.format(
            sum(zip(*number_of_samples_in_positive_and_negative_validation_sets)[0]) /
                float(len(number_of_samples_in_positive_and_negative_validation_sets)),
            sum(zip(*number_of_samples_in_positive_and_negative_validation_sets)[1]) /
                float(len(number_of_samples_in_positive_and_negative_validation_sets))))

        y_pred_all = np.array(y_pred_all, dtype=float)
        y_true_all = np.array(y_true_all, dtype=bool)
        dfs_left_out = pd.concat(dfs_left_out, ignore_index=True)

        if 'DrugPair' in self.predictor_df.columns:
            # Take the highest prediction for each drug pair
            temp = [ (list(self.predictor_df[self.predictor_df['DrugPair']==x[0]]['Type'])[0]=='Pos', max(x[1]),) for x in preds_per_drug_pair.items() ]
            temp = np.array(temp)
            y_true_all_2 = temp[:,0]
            y_pred_all_2 = temp[:,1]
        else:
            y_true_all_2 = None
            y_pred_all_2 = None

        # Fit the predictor using all availible data
        self.clf.fit(data, labels)
        return y_true_all, y_pred_all, y_true_all_2, y_pred_all_2, dfs_left_out


    def get_index_pos(self, target, target_pair):
        return target in target_pair


    def get_x_y(self, df, columns_to_keep):
        x = np.float64(df[columns_to_keep].values)
        y = df['Type'] == 'Pos'
        return x, y


    def _leave_drug_pair_out(self):
        """
        """
        pos_neg_indices = []
        # Get positive and negative indices for each drug pair
        if 'DrugPair' not in self.predictor_df.columns:
            logger.debug("Don't have a DrugPair column!'")
            logger.debug("Doing leave one target pair out xval...")

            for target_1, target_2 in self.predictor_df['TargetPair'].apply(lambda tp: tp.split('_')).values:
                index_pos = self.predictor_df['TargetPair'].apply(lambda tp: (target_1 in tp.split('_')) and (target_2 in tp.split('_')))
                index_neg = self.predictor_df['TargetPair'].apply(lambda tp: (target_1 not in tp.split('_')) and (target_2 not in tp.split('_')))
                pos_neg_indices.append([np.array(index_pos, dtype=bool), np.array(index_neg, dtype=bool)])

        else:
            logger.debug("Have a DrugPair column!")
            if 'TargetPair' in self.predictor_df.columns:
                logger.debug("Have a TargetPair column!")
                logger.debug("Doing leave one drug-pair/target-pair out xval...")
                have_target_pair = True
            else:
                logger.debug("Don't have a TargetPair column!")
                logger.debug("Doing leave one drug-pair out xval...")
                have_target_pair = False

            visited_drug_pairs = set()
            for drug_1, drug_2 in self.predictor_df['DrugPair'].apply(lambda dp: dp.split('_')).values:
                # Don't create a leave-out set for the save drug pair but different target pairs
                if (drug_1, drug_2) in visited_drug_pairs:
                    continue
                visited_drug_pairs.add((drug_1, drug_2,))

                index_pos = self.predictor_df['DrugPair'].apply(lambda dp: (drug_1 in dp.split('_')) and (drug_2 in dp.split('_'))).values
                drugs_to_leave_out = set([x for xx in self.predictor_df['DrugPair'].iloc[index_pos].apply(lambda dp: dp.split('_')).values for x in xx])

                if not have_target_pair:
                    index_neg = self.predictor_df['DrugPair'].apply(
                        lambda drug_pair:
                            (drug_pair.split('_')[0] not in drugs_to_leave_out) and
                            (drug_pair.split('_')[1] not in drugs_to_leave_out),
                        ).values
                else:
                    targets_to_leave_out = set([x for xx in self.predictor_df['TargetPair'].iloc[index_pos].apply(lambda tp: tp.split('_')).values for x in xx])
                    index_neg = self.predictor_df[['DrugPair', 'TargetPair']].apply(
                        lambda drug_pair, target_pair:
                            (drug_pair.split('_')[0] not in drugs_to_leave_out) and
                            (drug_pair.split('_')[1] not in drugs_to_leave_out) and
                            (target_pair.split('_')[0] not in targets_to_leave_out) and
                            (target_pair.split('_')[1] not in targets_to_leave_out),
                        axis=1).values
                pos_neg_indices.append([np.array(index_pos, dtype=bool), np.array(index_neg, dtype=bool)])

        return pos_neg_indices


    def join_drug_pair_indices(self, pos_neg_indices, n_folds):
        items_per_fold = int(math.ceil(float(len(pos_neg_indices))/n_folds))
        n_folds = int(math.ceil(float(len(pos_neg_indices))/items_per_fold))
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


###############################################################################
# Junk I don't really need anymore...


    def _join_groups_nfold(self, kf_xv, n_folds):
        kf_xv_2 = []
        kf = len(kf_xv) / n_folds + 1
        for i in range(0, len(kf_xv), kf):
            train_new = kf_xv[i][0]
            test_new = kf_xv[i][1]
            for j in range(1, kf):
                if i+j >= len(kf_xv):
                    continue
                train_new *= kf_xv[i+j][0]
                test_new += kf_xv[i+j][1]
            kf_xv_2.append([train_new, test_new])
        return kf_xv_2


    def leave_drug_target_out(self, df):
        set_of_drugs = set()
        __ = df['DrugPair'].apply( lambda x: set_of_drugs.update(x.split('_')) )
        pos_neg_indices = []
        for drug in set_of_drugs:
            index_pos = np.array(df['DrugPair'].apply(lambda x: drug in x), dtype=bool)
            index_neg = np.array(df['DrugPair'].apply(lambda x: not (drug in x)), dtype=bool)

            set_of_drugs_pos = set([drug])
            __ = df['DrugPair'].iloc[index_pos].apply(lambda x: set_of_drugs_pos.update(x.split('_')))
            for drug_pos in set_of_drugs_pos:
                index_neg *= np.array(df['DrugPair'].apply(lambda x: not (drug_pos in x)), dtype=bool)

            set_of_targets_pos = set()
            __ = df['TargetPair'].iloc[index_pos].apply(lambda x: set_of_targets_pos.update(x.split('_')))
            for target_pos in set_of_targets_pos:
                index_neg *= np.array(df['TargetPair'].apply(lambda x: not (target_pos in x)), dtype=bool)

            pos_neg_indices.append([index_neg, index_pos])
        return pos_neg_indices


    def leave_target_out(self, df):
        set_of_targets = set()
        __ = df['TargetPair'].apply(lambda x: set_of_targets.update(x.split('_')))
        pos_neg_indices = []
        for target in set_of_targets:
            index_pos = np.array(df['TargetPair'].apply(lambda x: target in x), dtype=bool)
            index_neg = np.array(df['TargetPair'].apply(lambda x: not (target in x)), dtype=bool)

            set_of_targets_pos = set([target])
            __ = df['TargetPair'].iloc[index_pos].apply( lambda x: set_of_targets_pos.update(x.split('_')) )
            for target_pos in set_of_targets_pos:
                index_neg *= np.array(df['TargetPair'].apply(lambda x: not (target_pos in x)), dtype=bool)

            pos_neg_indices.append([index_neg, index_pos])
        return pos_neg_indices


    def get_columns_to_keep(self, df, idxs):
        columns_to_keep = [df.column[idx] for idx in idxs]
        logger.debug(columns_to_keep)
        return columns_to_keep


    def discretize_feature(self, pd_series):
        quantiles = sp.percentile(pd_series, [100./6, 100./6*2, 100./6*3, 100./6*4, 100./6*5])
        def discretize(feature):
            for i in range(len(quantiles)):
                if feature < quantiles[i]:
                    discrete_feature = i
                    break
                else:
                    discrete_feature = len(quantiles)
            return discrete_feature
        return pd_series.apply(discretize)


    def export_data_for_cytoscape(self, predictor_name):
        predictor_df = self.predictor_df
        temp = predictor_df.DrugPair.apply(lambda x: x.split('_')[0] + '\tpp\t' + x.split('_')[1])
        temp = pd.concat([predictor_df.Type, temp], axis=1)
        temp.to_csv('/home/kimlab1/strokach/databases/chemical_interactions/version_1/drug_network_%s.csv' % predictor_name,
            sep='\t',
            header=False,
            index=True)
        temp = predictor_df.TargetPair.apply(lambda x: x.split('_')[0] + '\tpp\t' + x.split('_')[1])
        temp = pd.concat([predictor_df.Type, temp], axis=1)
        temp.to_csv('/home/kimlab1/strokach/databases/chemical_interactions/version_1/target_network_%s.csv' % predictor_name,
            sep='\t',
            header=False,
            index=True)

        predictor_df = predictor_2_df
        temp = predictor_df.DrugPair.apply(lambda x: x.split('_')[0] + '\tpp\t' + x.split('_')[1])
        temp = pd.concat([predictor_df.Type, temp], axis=1)
        temp.to_csv('/home/kimlab1/strokach/databases/chemical_interactions/version_1/drug_network_%s.csv' % predictor_name,
            sep='\t',
            header=False,
            index=True)
        temp = predictor_df.TargetPair.apply(lambda x: x.split('_')[0] + '\tpp\t' + x.split('_')[1])
        temp = pd.concat([predictor_df.Type, temp], axis=1)
        temp.to_csv('/home/kimlab1/strokach/databases/chemical_interactions/version_1/target_network_%s.csv' % predictor_name,
            sep='\t',
            header=False,
            index=True)
        #export_data_for_cytoscape(predictor_df, 'p1')
        #export_data_for_cytoscape(predictor_df, 'p2')

###############################################################################


#%%
if __name__ == '__main__':
    import argparse
    import sqlalchemy as sa
    from sklearn import metrics
    from sklearn.ensemble import GradientBoostingClassifier

    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_file', type=str)
    parser.add_argument('--path_to_data', nargs='?', type=str, default='/home/kimlab1/strokach/databases/chemical_interactions/')
    parser.add_argument('--clf_type', nargs='?', type=int, default=None)
    parser.add_argument('--n_folds', nargs='?', type=int, default=40)
    parser.add_argument('--loss', nargs='?', type=str, default='deviance')
    parser.add_argument('--n_estimators', type=int)
    parser.add_argument('--learning_rate', type=float)
    parser.add_argument('--max_depth', type=int)
    parser.add_argument('--subsample', type=float)
    parser.add_argument('--column_idxs', nargs='?', type=str, default='')
    parser.add_argument('--output_db', nargs='?', type=str, default='')
    args = parser.parse_args()


    # Prepare predictor parameters
    parameters = {
        'n_estimators': args.n_estimators,
        'learning_rate': args.learning_rate,
        'max_depth': args.max_depth,
        'subsample': args.subsample}

    additional_parameters = {
        'loss': args.loss,
        'random_state': random_state}

    clf_options = dict(parameters.items() + additional_parameters.items())

    logger.info('Initializing predictor...')
    clf = GradientBoostingClassifier(**clf_options)
    pred = Predictor(args.input_file, args.path_to_data, clf)


    # Prepare predictor data
    if args.column_idxs:
        """ Iterative feature addition
        """
        if args.clf_type == 1:
            columns_all = chemical_similarity_columns + chemical_other_columns + target_pairwise_columns
            columns_to_keep = [columns_all[int(i)] for i in args.column_idxs.split(',')]
        elif args.clf_type == 2 or args.clf_type == 3:
            columns_all = chemical_similarity_columns + target_pairwise_columns
            columns_to_keep = [columns_all[int(i)] for i in args.column_idxs.split(',')]
        elif args.clf_type == 4:
            columns_to_keep = chemical_similarity_columns
    else:
        """ Use all but the chemical similarity features
        """
        if args.clf_type == 0:
            columns_to_keep = []
        elif args.clf_type == 1:
            columns_to_keep = [] # chemical_other_columns + target_pairwise_columns
        elif args.clf_type == 2 or args.clf_type == 3:
            columns_to_keep = [] # target_pairwise_columns

    # Perform cross-validation
    logger.info('Cross-validating predictor...')
    y_true_all, y_pred_all, y_true_all_2, y_pred_all_2, dfs_left_out = pred.cross_validate_predictor(args.n_folds, columns_to_keep=columns_to_keep)


    # Calculate ROC-AUC, PR-AUC and AC using cross-validation results
    y_true_all, y_pred_all = y_true_all, y_pred_all
    roc_auc = metrics.roc_auc_score(y_true_all, y_pred_all)
    p, r, t = metrics.precision_recall_curve(y_true_all, y_pred_all)
    fpr, tpr, thr = metrics.roc_curve(y_true_all, y_pred_all)
    pr_auc = metrics.auc(r, p)
    f1_score = metrics.f1_score(y_true_all, y_pred_all>0.5)
    accuracy_score = metrics.accuracy_score(y_true_all, y_pred_all>0.5)
    print('%0.5f %0.5f %0.5f' % (roc_auc, pr_auc, accuracy_score))

    # Save results to a local sqlite database
    output_parameter = {
        'input_file': args.input_file,
        'clf_type': args.clf_type,
        'n_folds': args.n_folds,
        'roc_auc': roc_auc,
        'pr_auc': pr_auc,
        'accuracy_score': accuracy_score}
    gridsearch_params = pd.DataFrame(dict(parameters.items() + output_parameter.items()), index=[0])

    if args.output_db:
        output_db_filename = args.output_db
    else:
        output_db_filename =  path_to_data + 'version_2.1/gridsearch.db'.format(args.clf_type)
    engine = sa.create_engine('sqlite:///' + output_db_filename)
    gridsearch_params.to_sql('predictor_{}'.format(args.clf_type), engine, if_exists='append', index=False)



