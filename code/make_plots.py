# -*- coding: utf-8 -*-
"""
Created on Sun Dec 15 23:21:26 2013

@author: alexey
"""

import datetime
try:
    import cPickle as pickle
except ImportError:
    import pickle
import subprocess

import numpy as np
import scipy as sp
import pandas as pd
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

from sklearn import metrics, preprocessing, cross_validation
from sklearn.linear_model import LogisticRegression
from sklearn.externals import joblib
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.ensemble.partial_dependence import plot_partial_dependence

import chemical_interactions.code.chemical_interactions_v2 as ci
from chemical_interactions.code.class_ml import split_train_test
import chemical_interactions.code.constants as constants

pd.options.mode.chained_assignment = None
FILL_NA = True
random_state = np.random.RandomState(42)

version_suffix = '_v2'



#%%############################################################################
# Make ROC curves to compare different classifiers

n_folds = ci.n_folds
path_to_data = ci.path_to_data
input_files = ci.input_files
output_folder = ci.output_folder


#%%
#2014-04-17 best parameters

###### Predictor 1 ######s
clf_name = 'Predictor 1'
clf_type = 1
clf_idx = 0
plot_col_idx = 0
input_file = input_files[clf_idx]
clf_options = {
    'loss': 'deviance',
    'n_estimators': 900,
    'learning_rate': 0.01,
    'max_depth': 2,
    'subsample': 0.3,
    'random_state': random_state}
# NEW PARAMS
#clf_options = {
#    'loss': 'deviance',
#    'n_estimators': 100,
#    'learning_rate': 0.06,
#    'max_depth': 2,
#    'subsample': 0.7,
#    'random_state': random_state}


## p1hs
clf_name = 'Predictor 1 High Scored'
clf_type = 1
clf_idx = 0
plot_col_idx = 1
input_file = input_files[clf_idx+1]
clf_options = {
    'loss': 'deviance',
    'n_estimators': 700,
    'learning_rate': 0.05,
    'max_depth': 2,
    'subsample': 0.5,
    'random_state': random_state}



###### Predictor 2 ######
clf_name = 'Predictor 2'
clf_type = 2
clf_idx = 6
plot_col_idx = 2
input_file = input_files[clf_idx]
clf_options = {
    'loss': 'deviance',
    'n_estimators': 900,
    'learning_rate': 0.01,
    'max_depth': 4,
    'subsample': 0.7,
    'random_state': random_state}
## Try using clf params from p2hs
#clf_options = {
#    'loss': 'deviance',
#    'n_estimators': 700,
#    'learning_rate': 0.1,
#    'max_depth': 4,
#    'subsample': 0.7,
#    'random_state': random_state}

# p2hs
clf_name = 'Predictor 2 High Scored'
clf_type = 2
clf_idx = 6
plot_col_idx = 3
input_file = input_files[clf_idx+1]
clf_options = {
    'loss': 'deviance',
    'n_estimators': 700,
    'learning_rate': 0.1,
    'max_depth': 4,
    'subsample': 0.7,
    'random_state': random_state}



###### Predictor 0 ######
# NEW PARAMS
clf_name = 'Predictor 0'
clf_type = 0
clf_idx = 12
plot_col_idx = 4
input_file = input_files[clf_idx]
clf_options = {
    'loss': 'deviance',
    'n_estimators': 25,
    'learning_rate': 0.03,
    'max_depth': 2,
    'subsample': 0.3,
    'random_state': random_state}



#%%############################################################################
# Train the classifiers
#for n_folds in [10, 20, 30, 40, 60, 100]:
for n_folds in [60]:
    print('Initializing predictor...')
    clf = GradientBoostingClassifier(**clf_options)
    pred = ci.Predictor(input_file, path_to_data, clf)

    print('Cross-validating predictor using {}-fold xval...'.format(n_folds))
    y_true_all, y_pred_all, y_true_all_2, y_pred_all_2, dfs_left_out = pred.cross_validate_predictor(n_folds)
    predictor_info = (clf, y_true_all, y_pred_all, y_true_all_2, y_pred_all_2, dfs_left_out) + (clf_options,)

    output_filename = path_to_data + output_folder + clf_name.replace(' ', '') + '_new_params_{}'.format(n_folds) + '.pkl'
    pickle.dump(predictor_info, open(output_filename, 'wb'), pickle.HIGHEST_PROTOCOL)
    print("Saved output to: {}".format(output_filename))



#%% (Extracted from the predictor class) GARBAGE
# Compute ROC curve and AUC
mean_tpr = 0.0
mean_fpr = np.linspace(0, 1, 100)
fpr_all = []
tpr_all = []
roc_auc_all = []

fpr, tpr, thresholds = metrics.roc_curve(y_true_all, y_pred_all)
mean_tpr = sp.interp(mean_fpr, fpr, tpr)
mean_tpr[0] = 0.0
roc_auc = metrics.auc(fpr, tpr)

fpr_all.append(fpr)
tpr_all.append(tpr)
roc_auc_all.append(roc_auc)

print(y_true_all)
print(y_pred_all)

score_auc = metrics.roc_auc_score(y_true_all, y_pred_all)
score_fmeasure = metrics.f1_score(y_true_all, y_pred_all>0.5)
score_accuracy = np.sum( (y_pred_all>0.5) == y_true_all ) / float(len(y_true_all))
print('score_auc: %.4f' % score_auc)
print('score_fmeasure: %.4f' % score_fmeasure)
print('score_accuracy: %.4f' % score_accuracy)

# Calculate precision and recall with respect to each class
# AUC best way to evaluate performace on these
if False:
    print(mean_tpr)
    print(len(kf_xv))
    mean_tpr /= len(kf_xv)
    mean_tpr[-1] = 1.0
    mean_auc = metrics.auc(mean_fpr, mean_tpr)
    print("Mean ROC AUC: %0.4f" % (mean_auc*100))
    print("Mean ROC AUC: %0.4f" % (np.mean(score_auc)*100))
    print("Mean f1_score: %0.4f" % (np.mean(score_fmeasure)*100))
    print("Mean accuracy: %0.4f" % (np.mean(score_accuracy)*100))


predictor_0 = pd.read_csv(path_to_data + 'version_2.1/predictors_2/predictor_0.tsv', sep='\t')

predictor_0_high_scored = pd.read_csv(path_to_data + 'version_2.1/predictors_2/predictor_1_high_scored.tsv', sep='\t')[predictor_0.columns]
predictor_0_independent_validation = pd.read_csv(path_to_data + 'version_2.1/predictors_2/predictor_1_independent_validation.tsv', sep='\t')[predictor_0.columns]
predictor_0_all_unused = pd.read_csv(path_to_data + 'version_2.1/predictors_2/predictor_1_all_unused.tsv', sep='\t')[[c for c in predictor_0.columns if c != 'Type']]
predictor_0_all_unused_pairs = pd.read_csv(path_to_data + 'version_2.1/predictors_2/predictor_1_all_unused_pairs.tsv', sep='\t')[[c for c in predictor_0.columns if c != 'Type']]

predictor_0_high_scored.to_csv(path_to_data + 'version_2.1/predictors_2/predictor_0_high_scored.tsv', sep='\t', index=False)
predictor_0_independent_validation.to_csv(path_to_data + 'version_2.1/predictors_2/predictor_0_independent_validation.tsv', sep='\t', index=False)
predictor_0_all_unused.to_csv(path_to_data + 'version_2.1/predictors_2/predictor_0_all_unused.tsv', sep='\t', index=False)
predictor_0_all_unused_pairs.to_csv(path_to_data + 'version_2.1/predictors_2/predictor_0_all_unused_pairs.tsv', sep='\t', index=False)



#%%
engine = sa.create_engine('postgres://postgres:@192.168.6.19:5432/kimlab')
sql_query = """
select cid_1, cid_2, atc_similarity, morganfingerprintr2_tanimoto, side_effect_similarity
from chemical_interactions_v2.predictor_1_all_unused_pairs_3
join chemical_interactions_v2.drug_atc_similarity using (cid_1, cid_2)
join chemical_interactions_v2.drug_chemical_similarity using (cid_1, cid_2)
join chemical_interactions_v2.drug_side_effect_similarity using (cid_1, cid_2);
"""
all_drug_pairs = pd.read_sql_query(sql_query, engine)

all_drug_pairs_filename = 'version_2.1/predictors_2/all_drug_pairs.tsv'
all_drug_pairs_results_filename = 'version_2.1/reports/predictor_0_all_unused_pairs_3_scored.tsv'
all_drug_pairs_2 = all_drug_pairs.drop_duplicates()
all_drug_pairs_2['DrugPair'] = all_drug_pairs_2[['cid_1', 'cid_2']].apply(lambda x: 'CID{}_CID{}'.format(*x), axis=1)
all_drug_pairs_2.to_csv(path_to_data + all_drug_pairs_filename, sep='\t', index=False)
pred = ci.Predictor(all_drug_pairs_filename, path_to_data)
data_test, labels_test = get_data_and_labels(pred.predictor_df)
all_drug_pairs_2['probas_p0'] = clf.predict_proba(data_test)[:,1]
all_drug_pairs_2.to_csv(path_to_data + 'version_2.1/reports/predictor_0_all_unused_pairs_3_scored.tsv', index=False, sep='\t')



#%% Make plots of chemic
def make_plots_of_chemical_features(df_pos_confirmed, df_neg_confirmed):
    import seaborn as sns
    sns.set_context('talk', font_scale=1.2)
    
    path_to_output = '/home/kimlab1/strokach/working/chemical_interactions/results/14-11-07/'
    
    fg, ax = plt.subplots(figsize=(10,6))
    df_pos_confirmed['side_effect_similarity'].hist(range=(0,0.6), bins=10, ax=ax)
    df_neg_confirmed['side_effect_similarity'].hist(range=(0,0.6), bins=10, ax=ax, alpha=0.7)
    ax.set_xlabel('Side effect similarity')
    ax.set_ylabel('Number of drug pairs')
    ax.legend(['Confirmed positive', 'Confirmed negative'])
    plt.savefig(path_to_output + 'side_effect_similarity_hist.png', bbox_inches='tight', dpi=150)
    plt.savefig(path_to_output + 'side_effect_similarity.pdf', bbox_inches='tight')
    plt.savefig(path_to_output + 'side_effect_similarity.eps', bbox_inches='tight')
    
    
    fg, ax = plt.subplots(figsize=(10,6))
    df_pos_confirmed['chemical_similarity'].hist(range=(0,1), bins=10, ax=ax)
    df_neg_confirmed['chemical_similarity'].hist(range=(0,1), bins=10, ax=ax, alpha=0.7)
    ax.set_xlabel('Chemical similarity')
    ax.set_ylabel('Number of drug pairs')
    ax.legend(['Confirmed positive', 'Confirmed negative'])
    plt.savefig(path_to_output + 'chemical_similarity_hist.png', bbox_inches='tight', dpi=150)
    plt.savefig(path_to_output + 'chemical_similarity_hist.pdf', bbox_inches='tight')
    plt.savefig(path_to_output + 'chemical_similarity_hist.eps', bbox_inches='tight')
    
    
    fg, ax = plt.subplots(figsize=(10,6))
    df_pos_confirmed['atc_similarity'].hist(range=(0,5), bins=10, ax=ax)
    df_neg_confirmed['atc_similarity'].hist(range=(0,5), bins=10, ax=ax, alpha=0.7)
    ax.set_xlabel('ATC code similarity')
    ax.set_ylabel('Number of drug pairs')
    ax.legend(['Confirmed positive', 'Confirmed negative'])
    plt.savefig(path_to_output + 'atc_code_similarity_hist.png', bbox_inches='tight', dpi=150)
    plt.savefig(path_to_output + 'atc_code_similarity_hist.pdf', bbox_inches='tight')
    plt.savefig(path_to_output + 'atc_code_similarity_hist.eps', bbox_inches='tight')

    
    # Make a feature dependence plot
    features = [0, 2, 1]
    pred = ci.Predictor(input_file, path_to_data)
    data_train, labels_train = get_data_and_labels(pred.predictor_df)
    fg, ax = plt.subplots(figsize=(8,10), facecolor='white')
    plot_partial_dependence(
        clf, data_train, features, n_cols=2, percentiles=(0.01, 0.99),
        feature_names=['ATC similarity', 'Chemical similarity', 'Side effect similarity'],
        n_jobs=3, grid_resolution=100, ax=ax)
    plt.savefig(path_to_output + 'drug_pair_feature_importances.png', bbox_inches='tight', dpi=150)
    plt.savefig(path_to_output + 'drug_pair_feature_importances.pdf', bbox_inches='tight')
    plt.savefig(path_to_output + 'drug_pair_feature_importances.eps', bbox_inches='tight')


def score_predictions(y_true_all, y_pred_all):
    roc_auc = metrics.roc_auc_score(y_true_all, y_pred_all)
    p, r, t = metrics.precision_recall_curve(y_true_all, y_pred_all)
    fpr, tpr, thr = metrics.roc_curve(y_true_all, y_pred_all)
    pr_auc = metrics.auc(r, p)
    f1_score = metrics.f1_score(y_true_all, y_pred_all>0.5)
    accuracy_score = metrics.accuracy_score(y_true_all, y_pred_all>0.5)
    print '%0.5f %0.5f %0.5f' % (roc_auc, pr_auc, accuracy_score)
    
    
    

#%%############################################################################
# Load data from the trained classifiers

n_folds = 60
#output_filename = path_to_data + output_folder + clf_name.replace(' ', '') + '.pickle'
#output_filename = path_to_data + output_folder + clf_name.replace(' ', '') + 'leave_drug_out_xval.pickle'
#output_filename = path_to_data + output_folder + clf_name.replace(' ', '') + '_{}.pkl'.format(n_folds)
output_filename = path_to_data + output_folder + clf_name.replace(' ', '') + '_new_params_{}'.format(n_folds) + '.pkl'
with open(output_filename, 'rb') as fh:
    predictor_info = pickle.load(fh)
clf, y_true_all_1, y_pred_all_1, y_true_all_2, y_pred_all_2, dfs_left_out, args = predictor_info

# Calculate ROC-AUC, PR-AUC and AC using cross-validation results
score_predictions(y_true_all_1, y_pred_all_1)



#%%############################################################################
def get_confirmed_testfile(input_testfile, confirmed_drugpairs):
    input_file_data = pd.read_csv(path_to_data + input_testfile, sep='\t')
    input_file_data.columns = [constants.columnname_mapping.get(column, column) for column in input_file_data.columns]

    curated_drug_pairs = pd.read_csv(path_to_data + confirmed_drugpairs, sep='\t')
    curated_drug_pairs['cid_1'] = [int(drug_pair.split('_')[0][3:]) for drug_pair in curated_drug_pairs.DrugPair.values]
    curated_drug_pairs['cid_2'] = [int(drug_pair.split('_')[1][3:]) for drug_pair in curated_drug_pairs.DrugPair.values]
    assert sum([x[0] > x[1] for x in curated_drug_pairs[['cid_1', 'cid_2']].values]) == 0
    curated_drug_pairs = curated_drug_pairs[['cid_1','cid_2']]

    input_file_data = input_file_data.merge(curated_drug_pairs, on=['cid_1', 'cid_2'])
    confirmed_testfile = input_testfile + '.confirmed'
    input_file_data.to_csv(path_to_data + confirmed_testfile, sep='\t', index=False)
    return confirmed_testfile


def get_data_and_labels(predictor_df):
    columns_to_drop = ['Type', 'TargetPair', 'DrugPair']
    columns_to_keep = [col for col in predictor_df.columns if col not in columns_to_drop]
    data = np.float64(predictor_df[columns_to_keep].values)
    if 'Type' in predictor_df.columns:
        labels = predictor_df['Type'] == 'Pos'
    else:
        labels = np.zeros( (data.shape[0],1), dtype=bool )
    return data, labels


# Draw ROC curves with AUC values
current_date = datetime.datetime.now()
datestamp = '{:02d}-{:02d}-{:02d}/'.format(current_date.year, current_date.month, current_date.day)[2:]
path_to_output = '/home/kimlab1/strokach/working/chemical_interactions/results/' + datestamp
subprocess.check_call('mkdir -p {}'.format(path_to_output), shell=True)
def save_roc_as_csv(filename, tpr, fpr, thresholds):
    filename  ='_'.join(clf_name.split(' ')) + '-' + '_'.join(filename.split(' '))
    temp_df = pd.DataFrame({'true_positive_rate': tpr, 'false_positive_rate': fpr, 'thresholds':thresholds})
    temp_df.to_csv(path_to_output + filename + '.tsv', sep='\t', index=False)



#%%
#plt.figure(num=None, figsize=(8, 5), dpi=150, facecolor='w', edgecolor='k')
plt.figure(num=None, figsize=(8, 7))
take_max_over_target_pairs = False
plt.subplot(1, 1, 1)

label_string = '{}-fold cross-validation'.format(n_folds)
if not take_max_over_target_pairs:
    y_true_all, y_pred_all = y_true_all_1, y_pred_all_1
else:
    y_true_all, y_pred_all = y_true_all_2, y_pred_all_2
roc_auc = metrics.roc_auc_score(y_true_all, y_pred_all)
fpr, tpr, thresholds = metrics.roc_curve(y_true_all, y_pred_all)
plt.plot(fpr, tpr, 'r-',  lw=2, label='%s (area = %0.5f)' % (label_string, roc_auc)) #np.random.rand(3,1)
save_roc_as_csv(label_string, tpr, fpr, thresholds)
#plt.plot(fpr, tpr, 'b-',  label='%s (area = %0.5f)' % (clf_name, roc_auc), lw=2) #np.random.rand(3,1)
#plt.plot(mean_fpr, mean_tpr, 'b:',  label='%s old xval (area = %0.5f)' % (clf_name, score_auc), lw=2) #np.random.rand(3,1)


label_string = 'Positive test set vs. negative test set'
input_testfile = input_files[2+clf_idx]
pred = ci.Predictor(input_testfile, path_to_data)
df = pred.predictor_df.copy()
data_test, labels_test = get_data_and_labels(df)
df['probas'] = clf.predict_proba(data_test)[:,1]
df['labels'] = labels_test
if take_max_over_target_pairs:
    df = df.groupby('DrugPair').max()
roc_auc = metrics.roc_auc_score(df['labels'].values, df['probas'].values)
fpr, tpr, thresholds = metrics.roc_curve(df['labels'].values, df['probas'].values)
plt.plot(fpr, tpr, 'b-', lw=2, label='%s (area = %0.5f)' % (label_string, roc_auc))
df_pos = df[df['labels']==1]
df_neg = df[df['labels']==0]
save_roc_as_csv(label_string, tpr, fpr, thresholds)


label_string = 'Curated positive test set vs. negative test set'
confirmed_testfile = get_confirmed_testfile(
    input_testfile,
    'version_2.1/predictors_2/Curated_IndependentPositiveSet.tsv')
pred = ci.Predictor(confirmed_testfile, path_to_data)
df_pos_confirmed = pred.predictor_df.copy()
data_test, labels_test = get_data_and_labels(df_pos_confirmed)
df_pos_confirmed['probas'] = clf.predict_proba(data_test)[:,1]
df_pos_confirmed['labels'] = labels_test
if take_max_over_target_pairs:
    df_pos_confirmed = df_pos_confirmed.groupby('DrugPair').max()
df_merge = pd.concat([df_pos_confirmed, df_neg], ignore_index=True)
roc_auc = metrics.roc_auc_score(df_merge['labels'].values, df_merge['probas'].values)
fpr, tpr, thresholds = metrics.roc_curve(df_merge['labels'].values, df_merge['probas'].values)
plt.plot(fpr, tpr, 'y-', lw=2, label='%s (area = %0.5f)' % (label_string, roc_auc))
save_roc_as_csv(label_string, tpr, fpr, thresholds)


label_string = 'Curated positive test set vs. curated negative test set'
confirmed_testfile = get_confirmed_testfile(
    input_testfile,
    'version_2.1/predictors_2/Curated_IndependentNegativeSet.tsv')
pred = ci.Predictor(confirmed_testfile, path_to_data)
df_neg_confirmed = pred.predictor_df.copy()
data_test, labels_test = get_data_and_labels(df_neg_confirmed)
df_neg_confirmed['probas'] = clf.predict_proba(data_test)[:,1]
df_neg_confirmed['labels'] = labels_test
if take_max_over_target_pairs:
    df_neg_confirmed = df_neg_confirmed.groupby('DrugPair').max()
df_merge = pd.concat([df_pos_confirmed, df_neg_confirmed], ignore_index=True)
roc_auc = metrics.roc_auc_score(df_merge['labels'].values, df_merge['probas'].values)
fpr, tpr, thresholds = metrics.roc_curve(df_merge['labels'].values, df_merge['probas'].values)
plt.plot(fpr, tpr, 'm-', lw=2, label='%s (area = %0.5f)' % (label_string, roc_auc))
save_roc_as_csv(label_string, tpr, fpr, thresholds)


label_string = 'Positive test set vs. random pairs'
#input_file_test = input_files[4+clf_idx]
input_file_test = input_files[5+clf_idx]
pred = ci.Predictor(input_file_test, path_to_data)
df_neg_xval = pred.predictor_df.copy()
data_test, labels_test = get_data_and_labels(df_neg_xval)
df_neg_xval['probas'] = clf.predict_proba(data_test)[:,1]
df_neg_xval['labels'] = labels_test
if take_max_over_target_pairs:
    df_neg_xval = df_neg_xval.groupby('DrugPair').max()
df_merge = pd.concat([df_pos, df_neg_xval], ignore_index=True)
roc_auc = metrics.roc_auc_score(df_merge['labels'].values, df_merge['probas'].values)
fpr, tpr, thresholds = metrics.roc_curve(df_merge['labels'].values, df_merge['probas'].values)
plt.plot(fpr, tpr, 'g-', lw=2, label='%s (area = %0.5f)' % (label_string, roc_auc))
save_roc_as_csv(label_string, tpr, fpr, thresholds)


label_string = 'Confirmed positive test set vs. random pairs'
df_merge = pd.concat([df_pos_confirmed, df_neg_xval], ignore_index=True)
roc_auc = metrics.roc_auc_score(df_merge['labels'].values, df_merge['probas'].values)
fpr, tpr, thresholds = metrics.roc_curve(df_merge['labels'].values, df_merge['probas'].values)
plt.plot(fpr, tpr, 'c-', lw=2, label='%s (area = %0.5f)' % (label_string, roc_auc))
save_roc_as_csv(label_string, tpr, fpr, thresholds)
#make_plots_of_chemical_features(df_pos_confirmed, df_neg_xval)
score_predictions(df_merge['labels'].values, df_merge['probas'].values)


plt.legend(loc="lower right", prop={'size': 'medium'})
plt.xlabel('False positive rate', size='large')
plt.ylabel('True positive rate', size='large')
plt.title(clf_name + ('' if not take_max_over_target_pairs else ' best target pair'), size='x-large')

#plt.subplots_adjust(left=0.05, bottom=0.08, right=0.98, top=0.94, wspace=0.2, hspace=0.2)
plt.subplots_adjust(left=0.09, bottom=0.09, right=0.96, top=0.93, wspace=0.2, hspace=0.2)

plt.savefig(path_to_output + clf_name.lower().replace(' ', '_') + '.png', bbox_inches='tight', dpi=150)
plt.savefig(path_to_output + clf_name.lower().replace(' ', '_') + '.pdf', bbox_inches='tight')
plt.savefig(path_to_output + clf_name.lower().replace(' ', '_') + '.eps', bbox_inches='tight')




#%%############################################################################
# Original positive and negative training sets

input_testfile = input_files[2+clf_idx]
pred = ci.Predictor(input_testfile, path_to_data)
predictor_df_copy = pred.predictor_df.copy()
data_test, labels_test = get_data_and_labels(predictor_df_copy)
predictor_df_copy['probas'] = clf.predict_proba(data_test)[:,1]
predictor_df_copy['labels'] = labels_test

#proba_test_pos = predictor_df_copy[predictor_df_copy['labels']==1].groupby(['DrugPair']).max()['probas'].values
proba_test_pos = predictor_df_copy[predictor_df_copy['labels']==1]['probas'].values
plt.subplot(3, 4, 1 + plot_col_idx)
plt.hist(proba_test_pos, range=np.arange(-0.05, 1.05), bins=20, color='r', linewidth=2,
         label='Positive test set\n(n = {})'.format(proba_test_pos.shape[0]))
plt.xlim(-0.05,1.05)
plt.legend(loc='upper center', prop={'size':12})
plt.title(clf_name, size='x-large')

#proba_test_neg = predictor_df_copy[predictor_df_copy['labels']==0].groupby(['DrugPair']).max()['probas'].values
proba_test_neg = predictor_df_copy[predictor_df_copy['labels']==0]['probas'].values
plt.subplot(3, 4, 5 + plot_col_idx)
plt.hist(proba_test_neg, range=np.arange(-0.05, 1.05), bins=20, color='b', linewidth=2,
         label='Negative test set\n(n = {})'.format(proba_test_neg.shape[0]))
plt.xlim(-0.05,1.05)
plt.legend(loc='upper center', prop={'size':12})
plt.ylabel('Number of pairs', size='large')


# Confirmed positive training set
confirmed_testfile = get_confirmed_testfile(
    input_testfile,
    'version_2.1/predictors_2/Curated_IndependentPositiveSet.tsv')
pred = ci.Predictor(confirmed_testfile, path_to_data)
predictor_df_copy = pred.predictor_df.copy()
data_test, labels_test = get_data_and_labels(predictor_df_copy)
predictor_df_copy['probas'] = clf.predict_proba(data_test)[:,1]
predictor_df_copy['labels'] = labels_test

#proba_test_pos = predictor_df_copy[predictor_df_copy['labels']==1].groupby(['DrugPair']).max()['probas'].values
proba_test_pos = predictor_df_copy[predictor_df_copy['labels']==1]['probas'].values
plt.subplot(3, 4, 1 + plot_col_idx)
plt.hist(proba_test_pos, range=np.arange(-0.05, 1.05), bins=20, color='m', linewidth=2,
         label='Confirmed positive test set\n(n = {})'.format(proba_test_pos.shape[0]))
plt.xlim(-0.05,1.05)
plt.legend(loc='upper center', prop={'size':12})
plt.title(clf_name, size='x-large')

# Confirmed negative training set
confirmed_testfile = get_confirmed_testfile(
    input_testfile,
    'version_2.1/predictors_2/Curated_IndependentNegativeSet.tsv')
pred = ci.Predictor(confirmed_testfile, path_to_data)
predictor_df_copy = pred.predictor_df.copy()
data_test, labels_test = get_data_and_labels(predictor_df_copy)
predictor_df_copy['probas'] = clf.predict_proba(data_test)[:,1]
predictor_df_copy['labels'] = labels_test

#proba_test_pos = predictor_df_copy[predictor_df_copy['labels']==1].groupby(['DrugPair']).max()['probas'].values
proba_test_neg = predictor_df_copy[predictor_df_copy['labels']==0]['probas'].values
plt.subplot(3, 4, 5 + plot_col_idx)
plt.hist(proba_test_neg, range=np.arange(-0.05, 1.05), bins=20, color='c', linewidth=2,
         label='Confirmed positive test set\n(n = {})'.format(proba_test_neg.shape[0]))
plt.xlim(-0.05,1.05)
plt.legend(loc='upper center', prop={'size':12})
plt.title(clf_name, size='x-large')


# Randomly-chosen negative training set
input_file_test = input_files[4+clf_idx]
pred = ci.Predictor(input_file_test, path_to_data)
predictor_df_copy = pred.predictor_df.copy()
data_test, labels_test = get_data_and_labels(predictor_df_copy)
predictor_df_copy['probas'] = clf.predict_proba(data_test)[:,1]
predictor_df_copy['labels'] = labels_test

#proba_test_pos = predictor_df_copy.groupby(['DrugPair']).max()['probas'].values
proba_test = predictor_df_copy['probas'].values
plt.subplot(3, 4, 9 + plot_col_idx)
ph = plt.hist(proba_test, range=np.arange(-0.05, 1.05), bins=20, color='g', linewidth=2,
              label='All unused pairs in Stitch\n(n = {})'.format(proba_test.shape[0]))
plt.xlim(-0.05,1.05)
plt.legend(loc='upper center', prop={'size':12})
plt.xlabel('Probability of being synergistic', size='large')


###############################################################################

if True:
    ph = plt.hist(proba_test[:,1], range=np.arange(-0.05, 1.05), bins=50, normed=True, color=colors[2], linewidth=2)
else:
    x = np.linspace(-0.05, 1.05, 1000)
    y, bin_edges = np.histogram(proba_test[:,1], bins=100, density=True)
    #y = y / max(y)
    bincenters = 0.5*(bin_edges[1:]+bin_edges[:-1])
    f = interp1d(bincenters, y, kind='linear', fill_value=0, bounds_error=False)
    plt.plot(x, constants.smooth(f(x), smooth_factor, 'flat')[smooth_factor/2-1:-smooth_factor/2],
    color=colors[0], ls=lstyle, linewidth=2)


plt.legend(loc="lower right",prop={'size':12})
plt.title('Final Predictors 1 and 2')
plt.xlabel('False positive rate')
plt.ylabel('True positive rate')
plt.show()



def make_plot(feature_labels, feature_importances, title_string):
    fig = plt.figure(num=None, figsize=(16, 8), dpi=100, facecolor='w', edgecolor='k')
    bar_labels = feature_labels
    n_groups = len(feature_importances)
    # fig, ax = plt.subplots()
    index = np.arange(n_groups) + 0.35
    bar_width = 0.7
    opacity = 1 # 0.4
    # error_config = {'ecolor': '0.3'}
    # fig, ax = plt.subplots(figsize=(12,6), dpi=120)
    rects_1 = plt.bar(index, feature_importances, bar_width, alpha=opacity, color='b', figure=fig)

    plt.ylabel('Feature importance', size='large')
    plt.title(title_string, size='large')
    plt.xticks(index + bar_width/2, bar_labels, size='large', rotation=90, ha='center')
    # plt.xlim(0, 4.35)

    plt.tight_layout()
    plt.show()

ci.log.info(clf_type)
if clf_type == 1:
    columns_to_keep = ci.chemical_other_columns + ci.target_pairwise_columns
elif clf_type == 2 or clf_type == 3:
    columns_to_keep = ci.target_pairwise_columns

make_plot(columns_to_keep, clf.feature_importances_, clf_name)


###############################################################################
# Old code below:
# Load data from spreadsheet files
drug_data_pos = pd.read_csv('../drug_data_pos.tsv', sep='\t', index_col=0)
drug_data_neg = pd.read_csv('../drug_data_neg.tsv', sep='\t', index_col=0)
drug_data_columns = drug_data_pos.columns

target_data_pos = pd.read_csv('../target_data_pos.tsv', sep='\t', index_col=0)
target_data_neg = pd.read_csv('../target_data_neg.tsv', sep='\t', index_col=0)
target_data_columns = target_data_pos.columns

drug_data_merged_pos = pd.read_csv('../drug_data_merged_pos.tsv', sep='\t', index_col=0)
del drug_data_merged_pos['target_pair']
drug_data_merged_neg = pd.read_csv('../drug_data_merged_neg.tsv', sep='\t', index_col=0)
del drug_data_merged_neg['target_pair']
merged_data_columns = drug_data_merged_pos.columns

# Normalise columns to be between 0 and 1
mili_columns = ['functional_association',]
centi_columns = ['go_all_shared_fraction', 'go_bp_shared_fraction', 'go_cc_shared_fraction', 'go_mf_shared_fraction',]
for column in mili_columns:
    target_data_pos[column] = target_data_pos[column] / float(1000)
    target_data_neg[column] = target_data_neg[column] / float(1000)
    drug_data_merged_pos[column] = drug_data_merged_pos[column] / float(1000)
    drug_data_merged_neg[column] = drug_data_merged_neg[column] / float(1000)
for column in centi_columns:
    target_data_pos[column] = target_data_pos[column] / float(100)
    target_data_neg[column] = target_data_neg[column] / float(100)
    drug_data_merged_pos[column] = drug_data_merged_pos[column] / float(100)
    drug_data_merged_neg[column] = drug_data_merged_neg[column] / float(100)

# Deal with missing values
if FILL_NA:
    # Fill in missing chemical similarity features
    drug_data_pos.fillna(drug_data_pos.mean(), inplace=True)
    drug_data_neg.fillna(drug_data_neg.mean(), inplace=True)

    target_data_pos.fillna(target_data_pos.mean(), inplace=True)
    target_data_neg.fillna(target_data_neg.mean(), inplace=True)

    drug_data_merged_pos.fillna(drug_data_merged_pos.mean(), inplace=True)
    drug_data_merged_neg.fillna(drug_data_merged_neg.mean(), inplace=True)
else:
    # Remove rows with missing values
    drug_data_pos = drug_data_pos[~np.isnan(drug_data_pos).any(axis=1)]
    drug_data_neg = drug_data_neg[~np.isnan(drug_data_neg).any(axis=1)]

    target_data_pos = target_data_pos[~np.isnan(target_data_pos).any(axis=1)]
    target_data_neg = target_data_neg[~np.isnan(target_data_neg).any(axis=1)]

    drug_data_merged_pos = drug_data_merged_pos[~np.isnan(drug_data_merged_pos).any(axis=1)]
    drug_data_merged_neg = drug_data_merged_neg[~np.isnan(drug_data_merged_neg).any(axis=1)]

clf_drugs = GradientBoostingClassifier(loss='deviance',
                    n_estimators=2000,
                    learning_rate=0.08,
                    max_depth=4,
                    subsample=0.5,
                    random_state=random_state)
clf_targets = GradientBoostingClassifier(loss='deviance',
                    n_estimators=2000,
                    learning_rate=0.08,
                    max_depth=4,
                    subsample=0.5,
                    random_state=random_state)
clf_merged = GradientBoostingClassifier(loss='deviance',
                    n_estimators=2000,
                    learning_rate=0.08,
                    max_depth=4,
                    subsample=0.5,
                    random_state=random_state)

train_x, test_x, train_y, test_y = split_train_test(drug_data_pos, drug_data_neg, 16)
train_x_drugs = np.concatenate( (train_x, train_y), axis=0 )
train_y_drugs = np.concatenate( (np.ones(train_x.shape[0], dtype=int), np.zeros(train_y.shape[0], dtype=int)), axis=0 )
test_x_drugs = np.concatenate( (test_x, test_y), axis=0 )
test_y_drugs = np.concatenate( (np.ones(test_x.shape[0], dtype=int), np.zeros(test_y.shape[0], dtype=int)), axis=0 )
clf_drugs = clf_drugs.fit(train_x_drugs, train_y_drugs)

train_x, test_x, train_y, test_y = split_train_test(target_data_pos, target_data_neg, 10)
train_x_targets = np.concatenate( (train_x, train_y), axis=0 )
train_y_targets = np.concatenate( (np.ones(train_x.shape[0], dtype=int), np.zeros(train_y.shape[0], dtype=int)), axis=0 )
test_x_targets = np.concatenate( (test_x, test_y), axis=0 )
test_y_targets = np.concatenate( (np.ones(test_x.shape[0], dtype=int), np.zeros(test_y.shape[0], dtype=int)), axis=0 )
clf_targets = clf_targets.fit(train_x_targets, train_y_targets)

train_x, test_x, train_y, test_y = split_train_test(drug_data_merged_pos, drug_data_merged_neg, 16)
train_x_merged = np.concatenate( (train_x, train_y), axis=0 )
train_y_merged = np.concatenate( (np.ones(train_x.shape[0], dtype=int), np.zeros(train_y.shape[0], dtype=int)), axis=0 )
test_x_merged = np.concatenate( (test_x, test_y), axis=0 )
test_y_merged = np.concatenate( (np.ones(test_x.shape[0], dtype=int), np.zeros(test_y.shape[0], dtype=int)), axis=0 )
clf_merged = clf_merged.fit(train_x_merged, train_y_merged)


# Load classifiers

###############################################################################

# Drugs
probas_y_drugs = clf_drugs.predict_proba(test_x_drugs)[:,1]
fpr, tpr, thresholds = metrics.roc_curve(test_y_drugs, probas_y_drugs)
roc_auc = metrics.auc(fpr, tpr)
print "Area under the ROC curve : %f" % roc_auc
plt.clf()
plt.plot(fpr, tpr, label='SGB-DT Drug Pairs (Area = 0.89)', linewidth=2)

# Targets
probas_y_targets = clf_targets.predict_proba(test_x_targets)[:,1]
fpr, tpr, thresholds = metrics.roc_curve(test_y_targets, probas_y_targets)
roc_auc = metrics.auc(fpr, tpr)
print "Area under the ROC curve : %f" % roc_auc
plt.plot(fpr, tpr, label='SGB-DT Target Pairs (Area = 0.98)', linewidth=2)

# Merged
probas_y_merged = clf_merged.predict_proba(test_x_merged)[:,1]
fpr, tpr, thresholds = metrics.roc_curve(test_y_merged, probas_y_merged)
roc_auc = metrics.auc(fpr, tpr)
print "Area under the ROC curve : %f" % roc_auc
plt.plot(fpr, tpr, label='SGB-DT Merged (Area = 0.99)', linewidth=2)

plt.plot([0, 1], [0, 1], 'k--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.0])
plt.xlabel('False Positive Rate', size='large')
plt.ylabel('True Positive Rate', size='large')
#plt.title('Receiver operating characteristic', size='large')
plt.legend(loc="lower right")

plt.savefig('../course-project/ROC.png', dpi=600)


###############################################################################

# Drugs
probas_y_drugs = clf_drugs.predict_proba(test_x_drugs)[:,1]
precision, recall, thresholds = metrics.precision_recall_curve(test_y_drugs, probas_y_drugs)
area = metrics.auc(recall, precision)
print("Area Under Curve: %0.2f" % area)
plt.clf()
plt.plot(recall, precision, label='SGB-DT Drug Pairs (Area = %0.2f)' % area, linewidth=2)

# Targets
probas_y_targets = clf_targets.predict_proba(test_x_targets)[:,1]
precision, recall, thresholds = metrics.precision_recall_curve(test_y_targets, probas_y_targets)
area = metrics.auc(recall, precision)
print("Area Under Curve: %0.2f" % area)
plt.plot(recall, precision, label='SGB-DT Target Pairs (Area = %0.2f)' % area, linewidth=2)

# Merged
probas_y_merged = clf_merged.predict_proba(test_x_merged)[:,1]
precision, recall, thresholds = metrics.precision_recall_curve(test_y_merged, probas_y_merged)
area = metrics.auc(recall, precision)
print("Area Under Curve: %0.2f" % area)
plt.plot(recall, precision, label='SGB-DT Merged (Area = %0.2f)' % area, linewidth=2)


plt.xlabel('Recall')
plt.ylabel('Precision')
plt.ylim([0.0, 1.05])
plt.xlim([0.0, 1.0])
plt.title('Precision-Recall example: AUC=%0.2f' % area)
plt.legend(loc="lower left")
plt.show()

plt.savefig('../course-project/PRC.png', dpi=600)





