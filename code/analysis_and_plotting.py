# -*- coding: utf-8 -*-

import pickle
from sklearn import metrics

#temp = leave_target_out(predictor_df_targets)
#predictor_1_info = cross_validate_predictor(predictor_2_df_targets, columns_to_keep)

positive_target_pairs = negative_target_pairs = []
__ = map(positive_target_pairs.extend, predictor_df[predictor_df['Type']=='Pos']['TargetPair'].apply(lambda x: x.split('_')))
__ = map(negative_target_pairs.extend, predictor_df[predictor_df['Type']!='Pos']['TargetPair'].apply(lambda x: x.split('_')))
positive_target_pairs = set(positive_target_pairs)
negative_target_pairs = set(negative_target_pairs)

number_of_positives = sum(predictor_df['Type']=='Pos')
number_of_negatives = sum(predictor_df['Type']!='Pos')



###############################################################################
# Save data to get network classes




###############################################################################
path_to_data = '/home/kimlab1/strokach/working/chemical_interactions/data/version_2/'
with open(path_to_data + 'grid_search/Predictor_I_all_features.n_folds_30.loss_deviance.n_estimators_2000.learning_rate_0.080.max_depth_4.subsample_0.500.pickle') as fh: 
    predictor_1_info = pickle.load(fh)
with open(path_to_data + 'pickled/Predictor_with_drug_features_II.20xval_clf_info.pickle') as fh: 
    predictor_2_info = pickle.load(fh)
with open(path_to_data + 'pickled/PredictorII.20xval_clf_info_targets_only.pickle') as fh: 
    predictor_3_info = pickle.load(fh)


path_to_data = '/home/kimlab1/strokach/working/chemical_interactions/data/version_2/'
with open(path_to_data + 'pickled/PredictorHighScored_with_drug_features_I.20xval_clf_info.pickle') as fh: 
    predictor_1_info = pickle.load(fh)
with open(path_to_data + 'pickled/PredictorHighScored_with_drug_features_II.20xval_clf_info_targets_only.pickle') as fh: 
    predictor_2_info = pickle.load(fh)
with open(path_to_data + 'pickled/PredictorHighScored_II.20xval_clf_info_targets_only.pickle') as fh: 
    predictor_3_info = pickle.load(fh)


clf, y_true_all, y_pred_all, y_true_all_2, y_pred_all_2, __ = predictor_1_info

score_auc = metrics.roc_auc_score(y_true_all, y_pred_all)
score_fmeasure = metrics.f1_score(y_true_all, y_pred_all>0.5)
score_accuracy = np.sum( (y_pred_all>0.5) == y_true_all ) / float(len(y_true_all))
# Compute ROC curve and AUC
fpr, tpr, thresholds = metrics.roc_curve(y_true_all, y_pred_all)
mean_fpr = np.linspace(0, 1, 100)
mean_tpr = sp.interp(mean_fpr, fpr, tpr)
mean_tpr[0] = 0.0
roc_auc = metrics.auc(fpr, tpr)
plt.plot(mean_fpr, mean_tpr, 'b-',  label='Predictor 1 (area = %0.2f)' % score_auc, lw=2)

score_auc_2 = metrics.roc_auc_score(y_true_all_2, y_pred_all_2)
score_fmeasure_2 = metrics.f1_score(y_true_all_2, y_pred_all_2>0.5)
score_accuracy_2 = np.sum( (y_pred_all_2>0.5) == y_true_all_2 ) / float(len(y_true_all_2))
# Compute ROC curve and AUC
fpr, tpr, thresholds = metrics.roc_curve(y_true_all_2, y_pred_all_2)
mean_fpr = np.linspace(0, 1, 100)
mean_tpr = sp.interp(mean_fpr, fpr, tpr)
mean_tpr[0] = 0.0
roc_auc = metrics.auc(fpr, tpr)
plt.plot(mean_fpr, mean_tpr, 'r--',  label='Predictor 1 (area = %0.2f)' % score_auc_2, lw=2)

plt.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6), label='Random')         
plt.xlim([-0.05, 1.05])
plt.ylim([-0.05, 1.05])
plt.xlabel('False Positive Rate', fontsize=16)
plt.ylabel('True Positive Rate', fontsize=16)
plt.title(
    ' Predictor 1 individual: mean ROC: %.2f, mean f1 score: %.2f, mean accuracy: %.2f\n '
    % (score_auc, score_fmeasure, score_accuracy,) + 
    ' Predictor 1 maximum: mean ROC: %.2f, mean f1 score: %.2f, mean accuracy: %.2f\n '
    % (score_auc_2, score_fmeasure_2, score_accuracy_2,),
    fontsize=16)
plt.legend(loc="lower right")
plt.show()
  
        







if False:
    make_plot(
        columns_to_keep, 
        predictor_1_info[0].feature_importances_, 
        ' Predictor 1\n (mean ROC: %.2f, mean f1 score: %.2f, mean accuracy: %.2f)'
            % (predictor_1_info[3:6]))
    make_plot(
        columns_to_keep,
        predictor_2_info[0].feature_importances_, 
        ' Predictor 2 discretized\n (mean ROC: %.2f, mean f1 score: %.2f, mean accuracy: %.2f)'
            % (predictor_2_info[3:6]))
    
    plt.savefig('/home/kimlab1/strokach/documents/subgroup-meetings/drug-targets/140219/p1_5.png', dpi=300)


if False:
    fig, ax = plot_partial_dependence(predictor_2_info[0], 
        get_x_y(predictor_2_df_targets, columns_to_keep)[0], 
        features=range(12,23), 
        feature_names=columns_to_keep)




# Plot ROC curve
# Credit: http://scikit-learn.org/stable/auto_examples/plot_roc_crossval.html            
# for xv_idx, roc_auc in enumerate(roc_auc_all):
# plt.plot(fpr_all[xv_idx], tpr_all[xv_idx], lw=1, label='ROC fold %d (area = %0.02f)' % (xv_idx, roc_auc))
mean_fpr, mean_tpr, mean_auc = predictor_1_info[1], predictor_1_info[2], predictor_1_info[3]
plt.plot(mean_fpr, mean_tpr, 'b-',  label='Predictor 1 (area = %0.2f)' % mean_auc, lw=2)
#mean_fpr, mean_tpr, mean_auc = predictor_2_info[1], predictor_2_info[2], predictor_2_info[3]
#plt.plot(mean_fpr, mean_tpr, 'r--', label='Predictor 2 (area = %0.2f)' % mean_auc, lw=2)
#mean_fpr, mean_tpr, mean_auc = predictor_3_info[1], predictor_3_info[2], predictor_3_info[3]
#plt.plot(mean_fpr, mean_tpr, 'g.-', label='Predictor 3 (area = %0.2f)' % mean_auc, lw=2)    
plt.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6), label='Random')         
plt.xlim([-0.05, 1.05])
plt.ylim([-0.05, 1.05])
plt.xlabel('False Positive Rate', fontsize=16)
plt.ylabel('True Positive Rate', fontsize=16)
plt.title(
    ' Best targets\n '
    ' Predictor 1: mean ROC: %.2f, mean f1 score: %.2f, mean accuracy: %.2f\n '
    % (predictor_1_info[3:6]) + 
    ' Predictor 2: mean ROC: %.2f, mean f1 score: %.2f, mean accuracy: %.2f\n '
    % (predictor_2_info[3:6]) +
    ' Predictor 3: mean ROC: %.2f, mean f1 score: %.2f, mean accuracy: %.2f '
    % (predictor_3_info[3:6]),       
    fontsize=16)
plt.legend(loc="lower right")
plt.show()


# Plot colour diagrams of all features.
if False:
    data_pos = np.float64(predictor_2_df_targets[predictor_2_df_targets['Type']=='Pos'].values[:,2:])
    data_neg = np.float64(predictor_2_df_targets[predictor_2_df_targets['Type']!='Pos'].values[:1222,2:])
    data = np.vstack((data_pos,data_neg))
    
    max_min_scaler = preprocessing.MinMaxScaler()
    data_2 = max_min_scaler.fit_transform(data)
    
    h = plt.imshow(data_2, interpolation='none', cmap=plt.cm.spectral, aspect='auto')
    ax = plt.gca()
    ax.set_xticks(range(23))
    ax.set_xticklabels(predictor_2_df_targets.columns[2:], rotation=90)
    plt.grid(True)
    plt.colorbar()
    
    predictor_df_targets.columns[2:]

#    plt.savefig('/home/alexey/documents/presentations/subgroup-meetings/splicing/131204/%s-by-%s.png' % (mutation_type, binning_by), format='png', dpi=600)
#    plt.close()
