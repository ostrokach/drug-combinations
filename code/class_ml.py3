# -*- coding: utf-8 -*-

import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

from sklearn import utils
from sklearn import preprocessing
from sklearn import metrics
from sklearn import cross_validation

from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.grid_search import GridSearchCV

from sklearn import datasets
from scipy.stats.stats import pearsonr
from sklearn.metrics import mean_squared_error

from sklearn.externals import joblib

import pickle as pickle
random_state = np.random.RandomState(131213)

import os

###############################################################################

FILL_NA = True

NEG_SIZE_TRAIN = 16
NEG_SIZE_TEST = 10
SPLIT_DATA_BY = ['', 'cut', 'reshape', 'resample', 'cut1']
SPLIT_DATA_BY = SPLIT_DATA_BY[4]
O_NEG_RESAMPLED = 10


###############################################################################


def split_train_test(pos_data, neg_data, NEG_SIZE_TRAIN=NEG_SIZE_TRAIN):
    """ Split the data into a test set and a validation set
    """
    m_total_pos = pos_data.shape[0]
    m_total_neg = neg_data.shape[0]
    
    pos_data_train, pos_data_test = cross_validation.train_test_split(pos_data, test_size=0.20, random_state=random_state)
    del pos_data # don't have access to it anymore, agh!
    
    m_neg_train = pos_data_train.shape[0] * NEG_SIZE_TRAIN
    m_neg_test = pos_data_test.shape[0] * NEG_SIZE_TEST
    assert neg_data.shape[0] >= (m_neg_train + m_neg_test)
        
    # Split the negative data into training and validation sets
    neg_data = np.array(neg_data)
    neg_data = utils.shuffle(neg_data, random_state=random_state)
    
    neg_data_train = neg_data[m_neg_test:]
    neg_data_test = neg_data[:m_neg_test]   
    del neg_data # don't have access to it anymore, agh!
    
    o_neg_train = neg_data_train.shape[0] / m_neg_train
    assert neg_data_train.shape[0] >= o_neg_train*m_neg_train
    
    # Cut the negative training examples to be an exact multiple of the positive training examples
    if SPLIT_DATA_BY == 'cut' or SPLIT_DATA_BY == 'reshape':
        neg_data_train = neg_data_train[:o_neg_train*m_neg_train]
        
    # Split negative examples (for the training set) into o_neg_train-sized batches
    if SPLIT_DATA_BY == 'reshape' and o_neg_train >= 2:
        neg_data_train_new = np.empty( (m_neg_train,neg_data_train.shape[1],o_neg_train), dtype=float )
        for o_idx in range(o_neg_train):
            neg_data_train_new[:,:,o_idx] = neg_data_train[m_neg_train*o_idx:m_neg_train*(o_idx+1),:]
        neg_data_train = neg_data_train_new    
        assert neg_data_train.shape[0] == m_neg_train
    
    # Sample negative data to generate different training examples  
    if SPLIT_DATA_BY == 'resample':
        o_neg_train = O_NEG_RESAMPLED
        neg_data_train_new = np.empty( (m_neg_train,neg_data_train.shape[1],o_neg_train), dtype=float )
        for o_idx in range(o_neg_train):
            neg_data_train_new[:,:,o_idx] = utils.resample(neg_data_train, replace=True, n_samples=m_neg_train, random_state=random_state)
        neg_data_train = neg_data_train_new
        assert neg_data_train.shape[0] == m_neg_train
        
    print("Number of negative training datasets: %i" % o_neg_train)
            
    if SPLIT_DATA_BY == 'cut1':
        neg_data_train = neg_data_train[:m_neg_train]
    
    print("m_pos_train: %i, m_neg_train: %i, m_pos_total: %i, m_neg_total: %i" % \
    (pos_data_train.shape[0], neg_data_train.shape[0], m_total_pos, m_total_neg))
        
    print(pos_data_train.shape, pos_data_test.shape, neg_data_train.shape, neg_data_test.shape)
    
    return pos_data_train, pos_data_test, neg_data_train, neg_data_test


def run_ml_classifiers(data_set):
        
        x_train = data_set[0][0]
        y_train = data_set[0][1]
        
        classifiers = [
            ['LR (l1)', 
             [{'C': [80000,40000,20000,10000,5000,2000,1000,500,300,100,50,10,5,2,1,0.9,0.7,0.5,0.3,0.1,0.08,0.05,0.01],
              'class_weight':['auto', None],
              'penalty':['l1','l2']}],
            LogisticRegression(C=500, 
                               penalty='l1',
                               fit_intercept=True,
                               class_weight=None,
                               random_state=random_state)],
                               
#            ['LR (l2)',
#             [{'C': [40000,20000,10000,5000,4000,2000,1000,100,10,5,2,1,0.9,0.7,0.5,0.3,0.1,0.01],
#              'class_weight':['auto', None]}],
#            LogisticRegression(C=0.05, 
#                               penalty='l2',
#                               fit_intercept=True,
#                               class_weight=None)],
                               
            ['RF',
             [{'criterion': ['gini'],
               'n_estimators': [50,100,200,300,500,1000,2000,4000],
               'max_depth':[3,4,5,6,7,8,None],
               'max_features':['sqrt'],
               'min_samples_split':[1,2,3,4]}],
            RandomForestClassifier(n_estimators=500,
                               criterion='gini',
                               max_depth=None,
                               min_samples_split=2,
                               min_samples_leaf=1,
                               max_features='sqrt',
                               bootstrap=True,
                               oob_score=True,
                               random_state=random_state,
                               n_jobs=1,
                               compute_importances=None)],
                               
            ['GB',
            [{'n_estimators':[300,500,1000,2000,4000,6000],
             'learning_rate':[0.04,0.06,0.08],
             'max_depth':[2,3,4,6,None],
             'subsample':[0.5]}],
            GradientBoostingClassifier(loss='deviance',
                                n_estimators=2000, 
                                learning_rate=0.08, 
                                max_depth=4,
                                subsample=0.5,
                                random_state=random_state)]]
        
        # Type of cross-validation to perform
        loo_xv = cross_validation.LeaveOneOut(x_train.shape[0])
        kf_xv = cross_validation.StratifiedKFold(y_train, n_folds=6, indices=False)
     

#        for clf_idx, (clf_name, tuned_parameters, classifier) in enumerate(classifiers):
#            scores = ['f1']
#            for score in scores:
#                print "# Tuning %s hyper-parameters for %s" % (clf_name, score)
#                clf = GridSearchCV(classifier, tuned_parameters, cv=kf_xv, scoring=score, n_jobs=-1)
#                clf.fit(x_train, y_train)
#                print "Best parameters set found on development set:"
#                print clf.best_estimator_
#                print "Grid scores on development set:"
#                for params, mean_score, scores in clf.grid_scores_:
#                    print "%0.3f (+/-%0.03f) for %r" % (mean_score, scores.std() / 2, params)
#                print "Detailed classification report:" 
#                print "The model is trained on the full development set."
#                print "The scores are computed on the full evaluation set."
#                y_true, y_pred = data_set[1][1], clf.predict(data_set[1][0])
#                print metrics.classification_report(y_true, y_pred)
#                
#                # Save the best classifier
#                classifiers[clf_idx][2] = clf.best_estimator_


        for clf_idx, (clf_name, _, clf) in classifiers:
            # Compute score
            score_accuracy = []        
            score_auc = []
            score_fmeasure = []

            # Compute the AUC
            mean_tpr = 0.0
            mean_fpr = np.linspace(0, 1, 100)
            
            fpr_all = []
            tpr_all = []
            roc_auc_all = []
            y_pred_all = []
            for xv_idx, (train, test) in enumerate(kf_xv):
                
                clf.fit(x_train[train], y_train[train])
                score_accuracy.append(clf.score(x_train[test], y_train[test]))
                
                probas_ = clf.predict_proba(x_train[test])
                y_pred_all.append(probas_)

                score_auc.append(metrics.roc_auc_score(y_train[test], probas_[:,1]))
                score_fmeasure.append(metrics.f1_score(y_train[test], probas_[:,1]>0.5))
                
                # Compute ROC curve and AUC
                fpr, tpr, thresholds = metrics.roc_curve(y_train[test], probas_[:,1])
                mean_tpr += sp.interp(mean_fpr, fpr, tpr)
                mean_tpr[0] = 0.0
                roc_auc = metrics.auc(fpr, tpr)
                
                fpr_all.append(fpr)
                tpr_all.append(tpr)
                roc_auc_all.append(roc_auc)

            # Calculate precision and recall with respect to each class
            # AUC best way to evaluate performace on these
            mean_tpr /= len(kf_xv)
            mean_tpr[-1] = 1.0
            mean_auc = metrics.auc(mean_fpr, mean_tpr)
            print("Mean ROC AUC: %0.4f" % (mean_auc*100))
            print("Mean ROC AUC: %0.4f" % (np.mean(score_auc)*100))
            print("Mean f1_score: %0.4f" % (np.mean(score_fmeasure)*100))
            print("Mean accuracy: %0.4f" % (np.mean(score_accuracy)*100))
            # Plot ROC curve
            # Credit: http://scikit-learn.org/stable/auto_examples/plot_roc_crossval.html            
            if False:
                for xv_idx, roc_auc in enumerate(roc_auc_all):
                    plt.plot(fpr_all[xv_idx], tpr_all[xv_idx], lw=1, label='ROC fold %d (area = %0.02f)' % (xv_idx, roc_auc))
                plt.plot(mean_fpr, mean_tpr, 'k--', label='Mean ROC (area = %0.2f)' % mean_auc, lw=2)
                plt.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6), label='Luck')         
                plt.xlim([-0.05, 1.05])
                plt.ylim([-0.05, 1.05])
                plt.xlabel('False Positive Rate')
                plt.ylabel('True Positive Rate')
                plt.title('Receiver operating characteristic %s' % clf_name)
                plt.legend(loc="lower right")
                plt.show()
                input('Pause after plotting: ')

            
            ###################################################################
            
            
            clf.fit(x_train, y_train)
            
            print(clf_name + " training / cross-validation accuracy: %f / %f" % \
            (clf.score(x_train, y_train)*100,
            np.mean(cross_validation.cross_val_score(clf, x_train, y_train, cv=kf_xv, n_jobs=1))*100))
            
            if len(data_set) > 1:
                for idx, (x_data, y_data) in enumerate(data_set[1:]):
                    preds = clf.predict_proba(x_data)
                    print('*'*(idx+1) + clf_name +' ROC AUC on set ' + str(idx+1) + ': %f' % (metrics.roc_auc_score(y_data,preds[:,1])*100))
                    print('*'*(idx+1) + clf_name +' f1_score on set ' + str(idx+1) + ': %f' % (metrics.f1_score(y_data,preds[:,1]>0.5)*100))
                    print('*'*(idx+1) + clf_name +' accuracy on set ' + str(idx+1) + ': %f' % (metrics.accuracy_score(y_data,preds[:,1]>0.5)*100))
                    print('*'*(idx+1) + clf_name +' accuracy on set ' + str(idx+1) + ': %f' % (clf.score(x_data, y_data)*100))
                
            if clf_name[:2] == 'LR' and True:
                    hyperplane_coefs = clf.coef_[:]                    
                    pct_coefs_zero = np.mean(hyperplane_coefs == 0) * 100
                    print("Percentage of zero coefficients: %f" % pct_coefs_zero)


def run_ml_classifiers_final(train_x_all, train_y, test_x, test_y, exp_name):
        
        classifiers = [
            ['LR', 
             [],
            LogisticRegression(C=500, 
                               penalty='l1',
                               fit_intercept=True,
                               class_weight=None,
                               random_state=random_state)],
                               
            ['RF',
             [],
            RandomForestClassifier(n_estimators=500,
                               criterion='gini',
                               max_depth=None,
                               min_samples_split=2,
                               min_samples_leaf=1,
                               max_features='sqrt',
                               bootstrap=True,
                               oob_score=True,
                               random_state=random_state,
                               n_jobs=1,
                               compute_importances=None)],
                               
            ['GB',
            [],
            GradientBoostingClassifier(loss='deviance',
                                n_estimators=2000, 
                                learning_rate=0.08, 
                                max_depth=4,
                                subsample=0.5,
                                random_state=random_state)]]
                                
        clf_preds_all_overall = []
        for clf_name, _, clf in classifiers:
            print(clf_name + '-' + exp_name)
            
            if os.path.exists('../data/pickled/classifiers/' + clf_name + '-' + exp_name + '.pickle'):
                clfs_all = joblib.load('../data/pickled/classifiers/' + clf_name + '-' + exp_name + '.pickle')
            else:
                clfs_all = []
                for train_x in train_x_all:
                    clf.fit(train_x, train_y)
                    clfs_all.append(clf)
                joblib.dump(clfs_all, '../data/pickled/classifiers/' + clf_name + '-' + exp_name + '.pickle', compress=9)
            
            if len(clfs_all) == 1:
                    clf = clfs_all[0]
                    clf_preds = clf.predict_proba(test_x)[:,1]
                    clf_preds_all_overall.append(clf_preds)
                    print("AUC: %0.4f" % (metrics.roc_auc_score(test_y, clf_preds)*100))
                    print("f1_score on set: %0.6f" % (metrics.f1_score(test_y, clf_preds>0.5)*100))
                    print("accuracy on set: %0.6f" % (metrics.accuracy_score(test_y, clf_preds>0.5)*100))
            else:
                clf_preds_all = []
                for clf in clfs_all:
                    clf_preds_all.append(clf.predict_proba(test_x)[:,1])
                clf_preds = np.sum(clf_preds_all,axis=0) / len(clfs_all)
                clf_preds_all_overall.append(clf_preds)
                print("AUC: %0.4f" % (metrics.roc_auc_score(test_y, clf_preds)*100))
                print("f1_score on set: %0.6f" % (metrics.f1_score(test_y, clf_preds>0.5)*100))
                print("accuracy on set: %0.6f" % (metrics.accuracy_score(test_y, clf_preds>0.5)*100))  
        
        clf_preds_overall = np.sum(clf_preds_all_overall,axis=0) / len(classifiers)
        print("OVERALL CLF PRED RESULTS:")
        print("AUC: %0.4f" % (metrics.roc_auc_score(test_y, clf_preds_overall)*100))
        print("f1_score on set: %0.6f" % (metrics.f1_score(test_y, clf_preds_overall>0.5)*100))
        print("accuracy on set: %0.6f" % (metrics.accuracy_score(test_y, clf_preds_overall>0.5)*100))
        print("\n\n")

def run_ml(pos_data, neg_data, dataset_name):
    
    pos_data_train, pos_data_test, neg_data_train, neg_data_test = \
    split_train_test(pos_data, neg_data)
    
    test_x = np.concatenate( (pos_data_test, neg_data_test,), axis=0 )
    test_y = np.concatenate( (np.ones(pos_data_test.shape[0], dtype=int), 
                                   np.zeros(neg_data_test.shape[0], dtype=int)), 
                                   axis=0 )
    
    if len(neg_data_train.shape) == 2:
        # Only one negative training example
        train_x = np.concatenate( (pos_data_train, neg_data_train,), axis=0 )
        train_y = np.concatenate( (np.ones(pos_data_train.shape[0], dtype=int), 
                                   np.zeros(neg_data_train.shape[0], dtype=int)), 
                                   axis=0 )
#        run_ml_classifiers( [[train_x, train_y],] )
#        run_ml_classifiers( [[train_x, train_y],[test_x, test_y],] )
        if pos_data_train.shape[0] == neg_data_train.shape[0]:
            run_ml_classifiers_final([train_x], train_y, test_x, test_y, dataset_name + '-cut1')
        elif pos_data_train.shape[0] < neg_data_train.shape[0]:
            run_ml_classifiers_final([train_x], train_y, test_x, test_y, dataset_name + '-cut')
                                   
    
    else:
        # Multiple negative training examples, generated either by splitting or bootstrapping
        neg_list = np.array( list(range(np.size(neg_data_train,2))) )
        print(neg_list)
        train_x_all = []
        for neg_id in neg_list:
            print("Negative test set id: %i" % neg_id)
            
            train_x = np.concatenate( (pos_data_train, neg_data_train[:,:,neg_id],), axis=0 )
            train_y = np.concatenate( (np.ones(pos_data_train.shape[0], dtype=int),
                                       np.zeros(neg_data_train[:,:,neg_id].shape[0], dtype=int)), 
                                       axis=0 )
            train_x_all.append(train_x)
            
            non_train_id = neg_list[neg_list != neg_id]                       
            lo_x = neg_data_train[:,:,non_train_id[0]]
            for neg_id2 in non_train_id[1:]:
                lo_x = np.concatenate( (lo_x, neg_data_train[:,:,neg_id2]), axis=0 )
            lo_y = np.zeros(lo_x.shape[0], dtype=int)
            
#            run_ml_classifiers( [[train_x, train_y], [lo_x, lo_y],] )
#            run_ml_classifiers( [[train_x, train_y], [lo_x, lo_y], [test_x, test_y]] )         
        run_ml_classifiers_final(train_x_all, train_y, test_x, test_y, dataset_name + '-reshape')


if __name__ == '__main__':
       
    # Load data from spreadsheet files
    drug_data_pos = pd.read_csv('../drug_data_pos.tsv', sep='\t', index_col=0)
    drug_data_neg = pd.read_csv('../drug_data_neg.tsv', sep='\t', index_col=0)
    
    target_data_pos = pd.read_csv('../target_data_pos.tsv', sep='\t', index_col=0)
    target_data_neg = pd.read_csv('../target_data_neg.tsv', sep='\t', index_col=0)
    
    drug_data_merged_pos = pd.read_csv('../drug_data_merged_pos.tsv', sep='\t', index_col=0)
    del drug_data_merged_pos['target_pair']
    drug_data_merged_neg = pd.read_csv('../drug_data_merged_neg.tsv', sep='\t', index_col=0)
    del drug_data_merged_neg['target_pair']
    
    
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
    
    
    ###########################################################################
#    print "DRUG PAIRS!!!------------------------------------------------------"
#    run_ml(drug_data_pos, drug_data_neg, 'drug-pairs')
#    print "TARGET PAIRS!!!----------------------------------------------------"
#    run_ml(target_data_pos, target_data_neg, 'target-pairs')
    print("DRUG AND TARGET PAIRS MERGED!!!------------------------------------")
    run_ml(drug_data_merged_pos, drug_data_merged_neg, 'merged')



#scaler = preprocessing.StandardScaler(copy=False)
#scaler.fit(np.concatenate((drug_data_pos, drug_data_neg), axis=0))


#imp_pos = preprocessing.Imputer(missing_values=np.NaN, strategy='mean', axis=0)
#imp_pos.fit(chemical_similarity_pos_train)
#chemical_similarity_pos_train = imp_pos.transform(chemical_similarity_pos_train)
#chemical_similarity_pos_train = preprocessing.scale(chemical_similarity_pos_train)
#
#imp_neg = preprocessing.Imputer(missing_values=np.NaN, strategy='mean', axis=0)
#imp_neg.fit(chemical_similarity_neg_train)
#chemical_similarity_neg_train = imp_neg.transform(chemical_similarity_neg_train)
#chemical_similarity_neg_train = preprocessing.scale(chemical_similarity_neg_train)
#print np.mean(chemical_similarity_neg_train,0)
#
#chemical_similarity_pos_valid = preprocessing.scale(chemical_similarity_pos_valid)
#print np.mean(chemical_similarity_pos_valid,0)
#chemical_similarity_neg_valid = preprocessing.scale(chemical_similarity_neg_valid)
#print np.mean(chemical_similarity_neg_valid,0)


