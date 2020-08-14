# -*- coding: utf-8 -*-

import os
import time
import chemical_interactions_v2 as ci
import cPickle as pickle
import json
import subprocess
from sklearn import metrics

log = ci.log
header_columns = ci.header_columns
chemical_similarity_columns = ci.chemical_similarity_columns
chemical_other_columns = ci.chemical_other_columns
target_pairwise_columns = ci.target_pairwise_columns



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

##p2
#input_file = 'predictors/Predictor_with_drug_features_II.txt'
#clf_type = 2
#n_estimators = 500
#learning_rate = 0.05
#max_depth = 6
#subsample = 0.6
#output_folder = 'p2_iterative_feature_addition/'
#columns_all = chemical_similarity_columns + target_pairwise_columns


n_folds = 60

used_idxs = []
unused_idxs = range(len(columns_all))
scores = []
one_before = 0
two_before = 0
three_before = 0
#while not (three_before < two_before and two_before < one_before):
for i in range(len(columns_all)):
    idx_output_files = []
    output_file_scores = []
    for idx in unused_idxs:
        column_idxs = used_idxs + [idx]
        if 40 not in used_idxs:
            output_file = (path_to_data + output_folder + input_file.split('/')[-1].replace('.txt','') + 
            '.n_folds_%i.loss_%s.n_estimators_%i.learning_rate_%.3f.max_depth_%i.subsample_%.3f.column_idxs_%s.pickle' 
            % (n_folds, 'deviance', n_estimators, learning_rate, max_depth, subsample, ','.join([str(i) for i in column_idxs]),))
        else:
            output_file = (path_to_data + output_folder + input_file.split('/')[-1].replace('.txt','') + 
            '.n_folds_%i.loss_%s.n_estimators_%i.learning_rate_%.3f.max_depth_%i.subsample_%.3f.column_idxs_%s.pickle' 
            % (n_folds, 'deviance', n_estimators, learning_rate, max_depth, subsample, ','.join([str(i) for i in column_idxs[column_idxs.index(40):]]),))            
        idx_output_files.append([idx, output_file])
        if os.path.isfile(output_file):
            continue
        log.debug(output_file)
        system_command = (
            'submitjob 6 -m 3 -c 1 '
            'python2.7 chemical_interactions_v2.py '
            ' --path_to_data ' + path_to_data + 
            ' --input_file ' + input_file + 
            ' --output_file ' + output_file + 
            ' --n_estimators ' + str(n_estimators) + 
            ' --learning_rate ' + str(learning_rate) + 
            ' --max_depth ' + str(max_depth) + 
            ' --subsample ' + str(subsample) + 
            ' --clf_type ' + str(clf_type) + 
            ' --n_folds ' + str(n_folds) + 
            ' --loss deviance ' +
            ' --column_idxs ' + ','.join([str(i) for i in column_idxs]))
        log.debug(system_command)
        subprocess.check_call(system_command, shell=True)
    for idx, output_file in idx_output_files:
        while not os.path.isfile(output_file):
            log.debug('Waiting for file %s to finish...' % output_file)
            time.sleep(60)
        with open(output_file, 'rb') as fh:
            predictor_info = pickle.load(fh)
        __, y_true_all_1, y_pred_all_1, y_true_all_2, y_pred_all_2, args = predictor_info
        y_true_all, y_pred_all = y_true_all_2, y_pred_all_2
        score_auc = metrics.roc_auc_score(y_true_all, y_pred_all)
        score_fmeasure = metrics.f1_score(y_true_all, y_pred_all>0.5)
        score_accuracy = metrics.accuracy_score(y_true_all, y_pred_all>0.5)
        output_file_scores.append([idx, score_auc, score_fmeasure, score_accuracy])
    output_file_scores.sort(key=lambda x: x[1], reverse=True)
    used_idxs.append(output_file_scores[0][0])
    unused_idxs.remove(output_file_scores[0][0])
    scores.append(output_file_scores[0])
    if two_before:
        three_before = two_before
    if one_before:
        two_before = one_before
    one_before = output_file_scores[0][1]
        
log.debug('saving validation curve data...')
with open(path_to_data + output_folder + input_file.split('/')[-1].replace('.txt','.final_data.json'), 'w') as fh:
    json.dump(scores, fh)
log.debug('done!')
            


