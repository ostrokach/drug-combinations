# -*- coding: utf-8 -*-
"""
Created on Sun Dec 15 23:21:26 2013

@author: alexey
"""
from __future__ import print_function
import os
import numpy as np
import pandas as pd
import sqlalchemy as sa
import matplotlib.pyplot as plt
import cPickle as pickle
from sklearn import metrics


###################################################################################################
#%% Set parameters
try:
    code_path = os.path.dirname(os.path.realpath(__file__))
except NameError:
    code_path = '/home/kimlab1/strokach/working/chemical_interactions/code_release'

os.chdir(code_path)
import predictor

n_folds = 60
path_to_data = './'
input_folder = 'predictor_input/'
output_folder = 'predictor_output/'
input_files = [
    'predictor_1.tsv',
    'predictor_1_high_scored.tsv',
    'predictor_1_independent_validation.tsv',
    'predictor_1_all_unused.tsv',
    'predictor_1_all_unused_pairs.tsv',
    'predictor_2.tsv',
    'predictor_2_high_scored.tsv',
    'predictor_2_independent_validation.tsv',
    'predictor_2_all_unused.tsv',
    'predictor_2_all_unused_pairs.tsv'
]


#%% Select the predictor that you want to work with
#'predictor_1', 'predictor_1hs', 'predictor_2', 'predictor_2hs'
predictor_id = 'predictor_1'
predictor_parameters = predictor.predictor_parameters_all[predictor_id]


###################################################################################################
#%% Train the predictor
pred = predictor.Predictor(
    predictor.predictor_parameters_all[predictor_id]['input_file'],
    path_to_data + input_folder,
    path_to_data + output_folder)
pred.logger.info('Done initializing predictor!')

pred.logger.info('Cross-validating predictor...')
clf, y_true_all, y_pred_all, y_true_all_perdrugpair, y_pred_all_perdrugpair, dfs_left_out = \
    pred.cross_validate_predictor(predictor_parameters['clf_options'], n_folds)

#%%
output_filename = path_to_data + output_folder + predictor_id + '.pickle'
pred.logger.info('Saving cross-validation results to ' + output_filename)
predictor_info = (clf, y_true_all, y_pred_all, y_true_all_perdrugpair, y_pred_all_perdrugpair,) + (predictor_parameters,)
pickle.dump(predictor_info, open(output_filename, 'wb'), pickle.HIGHEST_PROTOCOL)


###################################################################################################
#%% Load data from the trained predictor
output_filename = path_to_data + output_folder + predictor_id + '.pickle'
predictor_info = pickle.load(open(output_filename))
clf, y_true_all, y_pred_all, y_true_all_perdrugpair, y_pred_all_perdrugpair, predictor_parameters = predictor_info


#%% Print basic statistics
roc_auc = metrics.roc_auc_score(y_true_all, y_pred_all)
p, r, t = metrics.precision_recall_curve(y_true_all, y_pred_all)
fpr, tpr, thr = metrics.roc_curve(y_true_all, y_pred_all)
pr_auc = metrics.auc(r, p)
f1_score = metrics.f1_score(y_true_all, y_pred_all>0.5)
accuracy_score = metrics.accuracy_score(y_true_all, y_pred_all>0.5)
print('ROC AUC: {:.5f}, PR AUC: {:.5f}, ACCURACY: {:.5f}'.format(roc_auc, pr_auc, accuracy_score))


#%% Use the trained classifier to score some test data
label_string = 'Positive test set vs. random pairs'
input_file_test = input_files[4+predictor_parameters['clf_idx']]
pred = predictor.Predictor(
    input_file_test,
    path_to_data + input_folder,
    path_to_data + output_folder)
data_test, labels_test = pred.get_data_and_labels()
df_neg_xval = pred.predictor_df.copy()
df_neg_xval['score_predictor_1'] = clf.predict_proba(data_test)[:,1]


#%% Load the trained high-scored classifier and run this code to get hs predictions
df_neg_xval['score_predictor_1_hs'] = clf.predict_proba(data_test)[:,1]


#%% Save results to a database
df_neg_xval = pd.read_csv(
    '/home/kimlab1/strokach/working/chemical_interactions/code_release/results/'
    'predictor_2_all_unused_pairs_scored.tsv', sep='\t')
del df_neg_xval['probas']
engine = sa.create_engine('postgresql://elaspic:elaspic@192.168.6.19:5432/kimlab')
meta = sa.MetaData(engine, schema='chemical_interactions_v2')
meta.reflect()
pdsql = pd.io.sql.PandasSQLAlchemy(engine, meta=meta)
pdsql.to_sql(df_neg_xval, 'predictor_2_all_unused_pairs_scored', if_exists='append', index=False)


#%%################################################################################################
# Draw ROC curves with AUC values for all classifiers
plt.figure(num=None, figsize=(8, 7))
take_max_over_target_pairs = False
plt.subplot(1, 1, 1)

label_string = '60-fold cross-validation'
if take_max_over_target_pairs:
    y_true_all, y_pred_all =  y_true_all_perdrugpair, y_pred_all_perdrugpair
roc_auc = metrics.roc_auc_score(y_true_all, y_pred_all)
fpr, tpr, thresholds = metrics.roc_curve(y_true_all, y_pred_all)
plt.plot(fpr, tpr, 'r-',  lw=2, label='%s (area = %0.5f)' % (label_string, roc_auc))


label_string = 'Positive test set vs. negative test set'
input_testfile = input_files[2+predictor_parameters['clf_idx']]
pred = predictor.Predictor(
    input_testfile,
    path_to_data + input_folder,
    path_to_data + output_folder)
df = pred.predictor_df.copy()
data_test, labels_test = pred.get_data_and_labels()
df['probas'] = clf.predict_proba(data_test)[:,1]
df['labels'] = labels_test
if take_max_over_target_pairs:
    df = df.groupby('DrugPair').max()
roc_auc = metrics.roc_auc_score(df['labels'].values, df['probas'].values)
fpr, tpr, thresholds = metrics.roc_curve(df['labels'].values, df['probas'].values)
plt.plot(fpr, tpr, 'b-', lw=2, label='%s (area = %0.5f)' % (label_string, roc_auc))
df_pos = df[df['labels']==1]
df_neg = df[df['labels']==0]


label_string = 'Curated positive test set vs. negative test set'
confirmed_testfile = predictor.get_confirmed_testfile(
    input_testfile,
    'Curated_IndependentPositiveSet.tsv',
    path_to_data + input_folder,
    '_positive')
pred = predictor.Predictor(
    confirmed_testfile,
    path_to_data + input_folder,
    path_to_data + output_folder)
df_pos_confirmed = pred.predictor_df.copy()
data_test, labels_test = pred.get_data_and_labels()
df_pos_confirmed['probas'] = clf.predict_proba(data_test)[:,1]
df_pos_confirmed['labels'] = labels_test
if take_max_over_target_pairs:
    df_pos_confirmed = df_pos_confirmed.groupby('DrugPair').max()
df_merge = pd.concat([df_pos_confirmed, df_neg], ignore_index=True)
roc_auc = metrics.roc_auc_score(df_merge['labels'].values, df_merge['probas'].values)
fpr, tpr, thresholds = metrics.roc_curve(df_merge['labels'].values, df_merge['probas'].values)
plt.plot(fpr, tpr, 'y-', lw=2, label='%s (area = %0.5f)' % (label_string, roc_auc))


label_string = 'Curated positive test set vs. curated negative test set'
confirmed_testfile = predictor.get_confirmed_testfile(
    input_testfile,
    'Curated_IndependentNegativeSet.tsv',
    path_to_data + input_folder,
    '_negative')
pred = predictor.Predictor(
    confirmed_testfile,
    path_to_data + input_folder,
    path_to_data + output_folder)
df_neg_confirmed = pred.predictor_df.copy()
data_test, labels_test = pred.get_data_and_labels()
df_neg_confirmed['probas'] = clf.predict_proba(data_test)[:,1]
df_neg_confirmed['labels'] = labels_test
if take_max_over_target_pairs:
    df_neg_confirmed = df_neg_confirmed.groupby('DrugPair').max()
df_merge = pd.concat([df_pos_confirmed, df_neg_confirmed], ignore_index=True)
roc_auc = metrics.roc_auc_score(df_merge['labels'].values, df_merge['probas'].values)
fpr, tpr, thresholds = metrics.roc_curve(df_merge['labels'].values, df_merge['probas'].values)
plt.plot(fpr, tpr, 'm-', lw=2, label='%s (area = %0.5f)' % (label_string, roc_auc))
print(label_string, df_pos_confirmed.shape, df_neg_confirmed.shape)


label_string = 'Positive test set vs. random pairs'
input_file_test = input_files[4+predictor_parameters['clf_idx']]
pred = predictor.Predictor(
    input_file_test,
    path_to_data + input_folder,
    path_to_data + output_folder)
df_neg_xval = pred.predictor_df.copy()
data_test, labels_test = pred.get_data_and_labels()
df_neg_xval['probas'] = clf.predict_proba(data_test)[:,1]
df_neg_xval['labels'] = labels_test
if take_max_over_target_pairs:
    df_neg_xval = df_neg_xval.groupby('DrugPair').max()
df_merge = pd.concat([df_pos, df_neg_xval], ignore_index=True)
roc_auc = metrics.roc_auc_score(df_merge['labels'].values, df_merge['probas'].values)
fpr, tpr, thresholds = metrics.roc_curve(df_merge['labels'].values, df_merge['probas'].values)
plt.plot(fpr, tpr, 'g-', lw=2, label='%s (area = %0.5f)' % (label_string, roc_auc))
print(label_string, df_pos.shape, df_neg_xval.shape)


label_string = 'Confirmed positive test set vs. random pairs'
df_merge = pd.concat([df_pos_confirmed, df_neg_xval], ignore_index=True)
roc_auc = metrics.roc_auc_score(df_merge['labels'].values, df_merge['probas'].values)
fpr, tpr, thresholds = metrics.roc_curve(df_merge['labels'].values, df_merge['probas'].values)
plt.plot(fpr, tpr, 'c-', lw=2, label='%s (area = %0.5f)' % (label_string, roc_auc))
print(label_string, df_pos_confirmed.shape, df_neg_xval.shape)


plt.legend(loc="lower right", prop={'size': 'large'})
plt.xlabel('False positive rate')#, size='large')
plt.ylabel('True positive rate')#, size='large')
plt.title(predictor_id + ('' if not take_max_over_target_pairs else ' best target pair'), size='x-large')

plt.tight_layout()

path_to_output = '/home/kimlab1/strokach/working/chemical_interactions/code_release/results/'
plt.savefig(path_to_output + predictor_id + '.png')
plt.savefig(path_to_output + predictor_id + '.pdf')
plt.savefig(path_to_output + predictor_id + '.eps')


#%%
from common import constants

pred = predictor.Predictor(
    predictor.predictor_parameters_all[predictor_id]['input_file'],
    path_to_data + input_folder,
    path_to_data + output_folder)
df_original = pred.predictor_df.copy()

def print_df_stats(df, df_id):
    print(
        '', constants.underline(df_id.capitalize() + ' positive:'),
        'all: {}, drug_pairs: {}, target_pairs: {}'.format(
            df[df['Type'] == 'Pos'].shape[0],
            df[df['Type'] == 'Pos'].drop_duplicates('DrugPair').shape[0],
            df[df['Type'] == 'Pos'].drop_duplicates('TargetPair').shape[0]),
        sep='\n'
    )
    print(
        '', constants.underline(df_id.capitalize() + ' negative:'),
        'all: {}, drug_pairs: {}, target_pairs: {}'.format(
            df[df['Type'] == 'Neg'].shape[0],
            df[df['Type'] == 'Neg'].drop_duplicates('DrugPair').shape[0],
            df[df['Type'] == 'Neg'].drop_duplicates('TargetPair').shape[0]),
        sep='\n'
    )

print_df_stats(df_original, predictor_id.capitalize())


#%%

input_testfile = input_files[2+predictor_parameters['clf_idx']]
pred = predictor.Predictor(
    input_testfile,
    path_to_data + input_folder,
    path_to_data + output_folder)
df = pred.predictor_df.copy()
print(label_string, df[df['labels'] == 1].shape, df[df['labels'] == 0].shape)


confirmed_testfile = predictor.get_confirmed_testfile(
    input_testfile,
    'Curated_IndependentPositiveSet.tsv',
    path_to_data + input_folder,
    '_positive')
pred = predictor.Predictor(
    confirmed_testfile,
    path_to_data + input_folder,
    path_to_data + output_folder)
print(label_string, df_pos_confirmed.shape, df_neg.shape)


# Confirmed positive test set
df_pos_confirmed

# Random pairs
df_neg_xval



#%%############################################################################
# Draw histograms of prediction scores for all classifiers
plt.figure(num=None, figsize=(24, 12))

input_testfile = input_files[2 + predictor_parameters['clf_idx']]
pred = predictor.Predictor(
    input_testfile,
    path_to_data + input_folder,
    path_to_data + output_folder)
predictor_df_copy = pred.predictor_df.copy()
data_test, labels_test = pred.get_data_and_labels()
predictor_df_copy['probas'] = clf.predict_proba(data_test)[:,1]
predictor_df_copy['labels'] = labels_test

#proba_test_pos = predictor_df_copy[predictor_df_copy['labels']==1].groupby(['DrugPair']).max()['probas'].values
proba_test_pos = predictor_df_copy[predictor_df_copy['labels']==1]['probas'].values
plt.subplot(3, 4, 1 + predictor_parameters['plot_col_idx'])
plt.hist(proba_test_pos, range=np.arange(-0.05, 1.05), bins=20, color='r', linewidth=2, alpha=0.8,
         label='Positive test set\n(n = {})'.format(proba_test_pos.shape[0]))
plt.xlim(-0.05,1)
plt.legend(loc='upper center', prop={'size':12})
plt.title(predictor_parameters['clf_name'], size='x-large')

#proba_test_neg = predictor_df_copy[predictor_df_copy['labels']==0].groupby(['DrugPair']).max()['probas'].values
proba_test_neg = predictor_df_copy[predictor_df_copy['labels']==0]['probas'].values
plt.subplot(3, 4, 5 + predictor_parameters['plot_col_idx'])
plt.hist(proba_test_neg, range=np.arange(-0.05, 1.05), bins=20, color='b', linewidth=2, alpha=0.8,
         label='Negative test set\n(n = {})'.format(proba_test_neg.shape[0]))
plt.xlim(-0.05,1)
plt.legend(loc='upper center', prop={'size':12})
plt.ylabel('Number of pairs', size='large')


# Confirmed positive training set
confirmed_testfile = predictor.get_confirmed_testfile(
    input_testfile,
    'Curated_IndependentPositiveSet.tsv',
    path_to_data + input_folder,
    '_positive')
pred = predictor.Predictor(
    confirmed_testfile,
    path_to_data + input_folder,
    path_to_data + output_folder)
data_test, labels_test = pred.get_data_and_labels()
predictor_df_copy = pred.predictor_df.copy()
predictor_df_copy['probas'] = clf.predict_proba(data_test)[:,1]
predictor_df_copy['labels'] = labels_test

#proba_test_pos = predictor_df_copy[predictor_df_copy['labels']==1].groupby(['DrugPair']).max()['probas'].values
proba_test_pos = predictor_df_copy[predictor_df_copy['labels']==1]['probas'].values
plt.subplot(3, 4, 1 + predictor_parameters['plot_col_idx'])
plt.hist(proba_test_pos, range=np.arange(-0.05, 1.05), bins=20, color='m', linewidth=2, alpha=0.8,
         label='Confirmed positive test set\n(n = {})'.format(proba_test_pos.shape[0]))
plt.xlim(-0.05,1)
plt.legend(loc='upper center', prop={'size':12})

# Confirmed negative training set
confirmed_testfile = predictor.get_confirmed_testfile(
    input_testfile,
    'Curated_IndependentNegativeSet.tsv',
    path_to_data + input_folder,
    '_negative')
pred = predictor.Predictor(
    confirmed_testfile,
    path_to_data + input_folder,
    path_to_data + output_folder)
data_test, labels_test = pred.get_data_and_labels()
predictor_df_copy = pred.predictor_df.copy()
predictor_df_copy['probas'] = clf.predict_proba(data_test)[:,1]
predictor_df_copy['labels'] = labels_test

#proba_test_pos = predictor_df_copy[predictor_df_copy['labels']==1].groupby(['DrugPair']).max()['probas'].values
proba_test_neg = predictor_df_copy[predictor_df_copy['labels']==0]['probas'].values
plt.subplot(3, 4, 5 + predictor_parameters['plot_col_idx'])
plt.hist(proba_test_neg, range=np.arange(-0.05, 1.05), bins=20, color='c', linewidth=2, alpha=0.8,
         label='Confirmed negative test set\n(n = {})'.format(proba_test_neg.shape[0]))
plt.xlim(-0.05,1)
plt.legend(loc='upper center', prop={'size':12})
plt.title(predictor_parameters['clf_name'], size='x-large')


# Randomly-chosen negative training set
input_file_test = input_files[4 + predictor_parameters['clf_idx']]
pred = predictor.Predictor(
    input_file_test,
    path_to_data + input_folder,
    path_to_data + output_folder)
data_test, labels_test = pred.get_data_and_labels()
predictor_df_copy = pred.predictor_df.copy()
predictor_df_copy['probas'] = clf.predict_proba(data_test)[:,1]
predictor_df_copy['labels'] = labels_test

#proba_test_pos = predictor_df_copy.groupby(['DrugPair']).max()['probas'].values
proba_test = predictor_df_copy['probas'].values
plt.subplot(3, 4, 9 + predictor_parameters['plot_col_idx'])
ph = plt.hist(proba_test, range=np.arange(-0.05, 1.05), bins=20, color='g', linewidth=2, alpha=0.8,
              label='All unused pairs in Stitch\n(n = {})'.format(proba_test.shape[0]))
plt.xlim(-0.05,1)
plt.legend(loc='upper center', prop={'size':12})
plt.title(predictor_parameters['clf_name'], size='x-large')
plt.xlabel('Probability of being synergistic', size='large')


plt.tight_layout()