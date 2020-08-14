# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 12:27:21 2014

@author: alexey
"""

import os
import time
import subprocess
import argparse
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
import cPickle as pickle
from sklearn import metrics

import chemical_interactions.code.chemical_interactions_v2 as ci

random_state = np.random.RandomState(42)


#%% Parse arguments and set parameters
parser = argparse.ArgumentParser()
parser.add_argument('--mode', nargs='?', type=str)
parser.add_argument('--clf_type', nargs='?', type=int)
parser.add_argument('--reverse_order', nargs='?', type=bool)

args = parser.parse_args()

n_folds = ci.n_folds
path_to_data = ci.path_to_data
input_files = ci.input_files
output_folder = ci.output_folder

path_to_code = '/home/kimlab1/strokach/working/chemical_interactions/code'
os.chdir(path_to_code)


#%%
n_estimators_list = [25, 50, 100, 200, 300, 400, 500, 600, 750, 1000, 1250, 1500, 1750, 2000, 2500, 3000]
learning_rate_list = [0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.10, 0.12, 0.15, 0.20, 0.40, 0.50, 0.60, 0.80, 1.00]
max_depth_list = [2, 3, 4, 5, 6]
subsample_list = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

#n_estimators_list = [300, 500, 700, 900]
#learning_rate_list = [0.01, 0.05, 0.1, 0.2]
#max_depth_list = [2, 4, 6]
#subsample_list = [0.3, 0.5, 0.7]
#max_features = ['sqrt', None]

#output_folder = 'results/grid_search_p%i_no_drug_similarity/' % args.clf_type
#input_files = [input_files[args.clf_type - 1], input_files[args.clf_type + 1]]

input_files = [input_files[10]]

if args.reverse_order:
    input_files = input_files.reverse()
    n_estimators_list = n_estimators_list.reverse()
    learning_rate_list = learning_rate_list.reverse()
    max_depth_list = max_depth_list.reverse()
    subsample_list = subsample_list.reverse()

if args.mode and args.mode == 'run':
    counter = 0
    for input_file in input_files:
        for n_estimators in n_estimators_list:
            for learning_rate in learning_rate_list:
                for max_depth in max_depth_list:
                    for subsample in subsample_list:
#                        output_file_suffix = \
#                            '.n_folds_%i.n_estimators_%i.learning_rate_%.3f.max_depth_%i.subsample_%.3f.clf_type_%i.pickle' \
#                            % (n_folds, n_estimators, learning_rate, max_depth, subsample, args.clf_type,)
#                        output_file = (output_folder + input_file.split('/')[-1].replace('.txt','') + output_file_suffix)
#                        if not os.path.isfile(output_file):
                        system_command = (
#                            ' submitjob 6 -m 3 -c 1 python python' +
                            ' echo ' +
                            ' --path_to_data ' + path_to_data +
                            ' --input_file ' + input_file +
#                            ' --output_file ' + output_file +
                            ' --clf_type ' + '{:d}'.format(args.clf_type) +
                            ' --n_folds ' + '{:d}'.format(n_folds) +
                            ' --n_estimators ' + '{:d}'.format(n_estimators) +
                            ' --learning_rate ' + '{:f}'.format(learning_rate) +
                            ' --max_depth ' + '{:d}'.format(max_depth) +
                            ' --subsample ' + '{:f}'.format(subsample))
                        subprocess.check_call(system_command, shell=True)
#                            counter += 1
#                            if (counter % 1000) == 0:
#                                time.sleep(3600)
    quit()





if args.mode and args.mode == 'parse':
    predictor_data = []
    for input_file in input_files:
        for n_estimators in n_estimators_list:
            for learning_rate in learning_rate_list:
                for max_depth in max_depth_list:
                    for subsample in subsample_list:
                        output_file_suffix = \
                            '.n_folds_%i.n_estimators_%i.learning_rate_%.3f.max_depth_%i.subsample_%.3f.clf_type_%i.pickle' \
                            % (n_folds, n_estimators, learning_rate, max_depth, subsample, args.clf_type,)
                        output_file = (output_folder + input_file.split('/')[-1].replace('.txt','') + output_file_suffix)
                        with open(path_to_data + output_file) as fh:
                            predictor_info = pickle.load(fh)
                        __, y_true_all, y_pred_all, y_true_all_2, y_pred_all_2, args = predictor_info
                        print args
                        score_auc = metrics.roc_auc_score(y_true_all, y_pred_all)
                        score_fmeasure = metrics.f1_score(y_true_all, y_pred_all>0.5)
                        score_accuracy = metrics.accuracy_score(y_true_all, y_pred_all>0.5)
                        predictor_data.append(
                            [args.input_file.split('/')[-1], args.input_file, args.output_file, args.clf_type,
                            args.n_folds, args.n_estimators, args.learning_rate, args.max_depth, args.subsample,
                            score_auc, score_fmeasure, score_accuracy])
    with open(path_to_data + output_folder + 'parsed_predictor_%i_output.pickle' % args.clf_type, 'wb') as fh:
        pickle.dump(predictor_data, fh)
    column_names = [
        'filename', 'input_filename', 'output_filename',
        'clf_type', 'n_folds', 'n_estimators', 'learning_rate', 'max_depth', 'subsample',
        'score_auc', 'score_fmeasure', 'score_accuracy']
    predictor_data_df = pd.DataFrame(predictor_data, columns=column_names)
    predictor_data_df.to_csv(path_to_data + output_folder + 'parsed_predictor_%i_output.tsv' % args.clf_type, index=False, sep='\t')
    quit()


###############################################################################



idx_of_interest = 4
drop_centrality_columns = 1
validatation_metric = ['score_auc_1', 'score_auc_2'][0]
best_params, predictor_info = best_parameters[(idx_of_interest,validatation_metric,)]
clf, y_true_all_1, y_pred_all_1, y_true_all_2, y_pred_all_2, args = predictor_info

clf_options = {'loss': args.loss,
                'n_estimators': args.n_estimators,
                'learning_rate': args.learning_rate,
                'max_depth': args.max_depth,
                'subsample': args.subsample,
                'random_state': chemical_interactions_v2.random_state}
if drop_centrality_columns:
    additional_columns_to_drop = \
        ['Degree',
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
else:
    additional_columns_to_drop = []
predictor = chemical_interactions_v2.Predictor(args.input_file, args.clf_type, args.path_to_data)
predictor_info = predictor.cross_validate_predictor(args.n_folds, 'drugs', clf_options, additional_columns_to_drop) + (args,)

output_file = (args.path_to_data + 'grid_search_2/' + args.input_file.split('/')[-1].replace('.txt','') +
                '.n_folds_%i.loss_%s.n_estimators_%i.learning_rate_%.3f.max_depth_%i.subsample_%.3f_drop_centrality_%i.pickle'
                    % (args.n_folds, args.loss, args.n_estimators, args.learning_rate, args.max_depth, args.subsample, drop_centrality_columns) )
with open(output_file, 'wb') as fh:
    pickle.dump(predictor_info, fh)

temp_1 = np.load('/home/kimlab1/strokach/working/databases/chemical_interactions/version_2/predictors/Predictor_with_drug_features_II.txt.1fold_xval.npy')
temp_30 = np.load('/home/kimlab1/strokach/working/databases/chemical_interactions/version_2/predictors/Predictor_with_drug_features_II.txt.30fold_xval.npy')
temp_60 = np.load('/home/kimlab1/strokach/working/databases/chemical_interactions/version_2/predictors/Predictor_with_drug_features_II.txt.60fold_xval.npy')
temp_100 = np.load('/home/kimlab1/strokach/working/databases/chemical_interactions/version_2/predictors/Predictor_with_drug_features_II.txt.100fold_xval.npy')




###############################################################################

if False:
    """ Plot chemical similarity features only
    """
    input_file = 'version_2/predictors/Predictor_with_drug_features_II.txt'
    path_to_data = '/home/kimlab1/strokach/working/databases/chemical_interactions/'

    args.clf_type = 0

    loss = 'deviance'
    n_estimators = 500
    learning_rate = 0.05
    max_depth = 6
    subsample = 0.6

    #chem_sim_columns = [''] + [str(x) for x in range(23)]
    chem_sim_columns = ['']

    plt.figure(num=None, figsize=(8, 6), dpi=150, facecolor='w', edgecolor='k')

    n_folds = 60

    for chem_sim_column in chem_sim_columns:
        output_file_suffix = \
            '.n_folds_%i.loss_%s.n_estimators_%i.learning_rate_%.3f.max_depth_%i.subsample_%.3f.column_idxs_%s.clf_type_%i.pickle' \
            % (n_folds, loss, n_estimators, learning_rate, max_depth, subsample, chem_sim_column, args.clf_type,)
        # print output_file_suffix
        output_file = (path_to_data + 'results/' + input_file.split('/')[-1].replace('.txt','') + output_file_suffix)
        with open(output_file, 'rb') as fh:
            predictor_info = pickle.load(fh)
        clf, y_true_all_1, y_pred_all_1, y_true_all_2, y_pred_all_2, args = predictor_info
        y_true_all, y_pred_all = y_true_all_2, y_pred_all_2
        score_auc = metrics.roc_auc_score(y_true_all, y_pred_all)
        score_fmeasure = metrics.f1_score(y_true_all, y_pred_all>0.5)
        score_accuracy = metrics.accuracy_score(y_true_all, y_pred_all>0.5)
        print '%s %0.5f %0.5f %0.5f' % (chem_sim_column, score_auc, score_fmeasure, score_accuracy)
        # Compute ROC curve and AUC
        fpr, tpr, thresholds = metrics.roc_curve(y_true_all, y_pred_all)
        mean_fpr = np.linspace(0, 1, 100)
        mean_tpr = sp.interp(mean_fpr, fpr, tpr)
        mean_tpr[0] = 0.0
        plt.plot(mean_fpr, mean_tpr, 'g-',  label='SGB %i-fold xval %s (area = %0.5f)' % (n_folds, str(chem_sim_column), score_auc), lw=2) #np.random.rand(3,1)


    plt.legend(loc="lower right",prop={'size':12})
    plt.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6), label='Random')
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.title('Chemical similarity features', size='large')
    plt.xlabel('False Positive Rate', size='medium')
    plt.ylabel('True Positive Rate', size='medium')
    plt.show()

###############################################################################


###############################################################################

if False:
    """ Make bar graphs showing the results of iterative feature addition
    """
    import json
    import matplotlib.pyplot as plt
    import chemical_interactions_v2 as ci

    log = ci.log
    header_columns = ci.header_columns
    chemical_similarity_columns = ci.chemical_similarity_columns
    chemical_other_columns = ci.chemical_other_columns
    target_pairwise_columns = ci.target_pairwise_columns
    all_columns = chemical_similarity_columns + chemical_other_columns + target_pairwise_columns
    path_to_data = '/home/kimlab1/strokach/working/databases/chemical_interactions/version_2/'
    loss = 'deviance'

    #p1
    input_file = 'predictors/Predictor_with_drug_features_I.txt'
    clf_type = 1
    n_estimators = 500
    learning_rate = 0.06
    max_depth = 2
    subsample = 0.4
    output_folder = 'p1_iterative_feature_addition/'
    columns_all = chemical_similarity_columns + chemical_other_columns + target_pairwise_columns
    with open(path_to_data + output_folder + input_file.split('/')[-1].replace('.txt','.final_data.json')) as fh:
        p_1 = json.load(fh)
    p_1_idxs, p_1_auc = zip(*p_1)[:2]

    #p2
    input_file = 'predictors/Predictor_with_drug_features_II.txt'
    clf_types = 2
    n_estimators = 500
    learning_rate = 0.05
    max_depth = 6
    subsample = 0.6
    output_folder = 'p2_iterative_feature_addition/'
    columns_all = chemical_similarity_columns + target_pairwise_columns
    with open(path_to_data + output_folder + input_file.split('/')[-1].replace('.txt','.final_data.json')) as fh:
        p_2 = json.load(fh)
    p_2_idxs, p_2_auc = zip(*p_2)[:2]

    title_string = 'Predictor 1'
    p_auc = p_1_auc
    p_idxs = p_1_idxs
    columns_all = chemical_similarity_columns + chemical_other_columns + target_pairwise_columns
    line_type = 'b-'

    title_string = 'Predictor 2'
    p_auc = p_2_auc
    p_idxs = p_2_idxs
    columns_all = chemical_similarity_columns + target_pairwise_columns
    line_type = 'r-'

    convert_labels = {
        'EB_Minimum': 'EB_Minimum_PPI',
        'EB_Maximum': 'EB_Maximum_PPI',
        'EB_Mean': 'EB_Mean_PPI',
        'EB_Fraction': 'EB_Fraction_PPI',
        'ShortestPathLength': 'ShortestPathLength_STRINGTopo',
    }

    plt.figure(num=None, figsize=(12, 6), dpi=150, facecolor='w', edgecolor='k')
    plt.plot(range(len(p_auc)), [1-auc for auc in p_auc], line_type, linewidth=3)
    plt.title(title_string, size='x-large')
    plt.ylabel('1 - AUC', size='large')
    plt.xticks(range(len(p_auc)), [convert_labels.get(columns_all[idx],columns_all[idx]) for idx in p_idxs], size='small', rotation=90, ha='center')
    plt.xlim(0, len(p_auc))
    plt.tight_layout()

###############################################################################




features_df = pd.read_csv(path_to_data + input_files[idx_of_interest], sep='\t')
if drop_centrality_columns:
    columns = [ column for column in list(features_df.columns)[3:] if column not in additional_columns_to_drop]
else:
    columns = list(features_df.columns)[3:]
make_plot(columns, clf.feature_importances_, '')

predictor_names = [
'Predictor 1 all targets',
'Predictor 1 best targets',
'Predictor 2 all targets',
'Predictor 2 best targets',
'Predictor 3 all targets',
'Predictor 3 best targets']
c = ['b', 'g', 'r', 'c', 'm', 'y']

if validatation_metric == 'score_auc_1':
    y_true_all, y_pred_all = y_true_all_1, y_pred_all_1
    score_auc, score_fmeasure, score_accuracy = best_params.score_auc_1, best_params.score_fmeasure_1, best_params.score_accuracy_1
    line_style = c[i] + '-'
elif validatation_metric == 'score_auc_2':
    y_true_all, y_pred_all = y_true_all_2, y_pred_all_2
    score_auc, score_fmeasure, score_accuracy = best_params.score_auc_2, best_params.score_fmeasure_2, best_params.score_accuracy_2
    line_style = c[i] + '--'

# Compute ROC curve and AUC
fpr, tpr, thresholds = metrics.roc_curve(y_true_all, y_pred_all)
mean_fpr = np.linspace(0, 1, 100)
mean_tpr = sp.interp(mean_fpr, fpr, tpr)
mean_tpr[0] = 0.0
plt.plot(mean_fpr, mean_tpr, line_style,  label='%s (area = %0.2f)' % (predictor_names[idx_of_interest], score_auc), lw=2)

plt.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6), label='Random')
plt.xlim([-0.05, 1.05])
plt.ylim([-0.05, 1.05])
plt.xlabel('False Positive Rate', fontsize=16)
plt.ylabel('True Positive Rate', fontsize=16)
#plt.title(
#    ' Predictor 1 individual: mean ROC: %.2f, mean f1 score: %.2f, mean accuracy: %.2f\n '
#    % (score_auc, score_fmeasure, score_accuracy,) +
#    ' Predictor 1 maximum: mean ROC: %.2f, mean f1 score: %.2f, mean accuracy: %.2f\n '
#    % (score_auc_2, score_fmeasure_2, score_accuracy_2,),
#    fontsize=16)
plt.legend(loc="lower right")
plt.show()

predictor_names = [
'Predictor 1 all targets',
'Predictor 1 best targets',
'Predictor 2 all targets',
'Predictor 2 best targets',
'Predictor 3 all targets',
'Predictor 3 best targets']
c = ['b', 'g', 'r', 'c', 'm', 'y']
plt.figure(num=None, figsize=(8, 6), dpi=150, facecolor='w', edgecolor='k')
for i in range(4,5):
    for j in ['score_auc_1', 'score_auc_2'][:2]:
        best_params, predictor_info = best_parameters[(i,j,)]
        clf, y_true_all_1, y_pred_all_1, y_true_all_2, y_pred_all_2, __ = predictor_info

        if j == 'score_auc_1':
            y_true_all, y_pred_all = y_true_all_1, y_pred_all_1
            score_auc, score_fmeasure, score_accuracy = best_params.score_auc_1, best_params.score_fmeasure_1, best_params.score_accuracy_1
            line_style = c[i] + '-'
        elif j == 'score_auc_2':
            y_true_all, y_pred_all = y_true_all_2, y_pred_all_2
            score_auc, score_fmeasure, score_accuracy = best_params.score_auc_2, best_params.score_fmeasure_2, best_params.score_accuracy_2
            line_style = c[i] + '--'

        # Compute ROC curve and AUC
        fpr, tpr, thresholds = metrics.roc_curve(y_true_all, y_pred_all)
        mean_fpr = np.linspace(0, 1, 100)
        mean_tpr = sp.interp(mean_fpr, fpr, tpr)
        mean_tpr[0] = 0.0
        plt.plot(mean_fpr, mean_tpr, line_style,  label='%s (area = %0.2f)' % (predictor_names[i], score_auc), lw=2)

plt.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6), label='Random')
plt.xlim([-0.05, 1.05])
plt.ylim([-0.05, 1.05])
plt.xlabel('False Positive Rate', fontsize=16)
plt.ylabel('True Positive Rate', fontsize=16)
#plt.title(
#    ' Predictor 1 individual: mean ROC: %.2f, mean f1 score: %.2f, mean accuracy: %.2f\n '
#    % (score_auc, score_fmeasure, score_accuracy,) +
#    ' Predictor 1 maximum: mean ROC: %.2f, mean f1 score: %.2f, mean accuracy: %.2f\n '
#    % (score_auc_2, score_fmeasure_2, score_accuracy_2,),
#    fontsize=16)
plt.legend(loc="lower right")
plt.show()

###############################################################################



