# -*- coding: utf-8 -*-
import pandas as pd
import cPickle as pickle

import numpy as np
import matplotlib.pyplot as plt

#
#h = plt.figure()
#n, bins_pos = np.histogram(p_db.target_data_pos['coexpression'], density=True, bins=np.arange(-0.7,1.1,0.1))
#_, bins_neg = np.histogram(p_db.target_data_neg['coexpression'], density=True, bins=np.arange(-0.7,1.1,0.1))
#
#mu, sigma = 200, 25
#x = mu + sigma*P.randn(10000)

# the histogram of the data with histtype='step'
plt.figure()
n, bins, patches = plt.hist(p_db.target_data_pos['coexpression'], 
                          bins=40, 
                          range=[-0.65,1.05],
                          normed=1, 
                          histtype='stepfilled',
                          facecolor='r',
                          alpha=0.75)
     

                   
n, bins, patches = plt.hist(p_db.target_data_neg['coexpression'], 
                          bins=40, 
                          range=[-0.65,1.05],
                          normed=1, 
                          histtype='stepfilled',
                          facecolor='b',
                          alpha=0.75)
plt.xlabel('Coexpression (Pearson R)', size='large')
plt.ylabel('Number of target pairs (normalised)', size='large')
plt.show()
plt.legend(['positive pairs', 'negative pairs'])


coexpression_data = np.concatenate((p_db.target_data_pos['coexpression'], p_db.target_data_neg['coexpression']))
coexpression_data = coexpression_data[~np.isnan(coexpression_data)]
coexpression_data_pos = np.array(p_db.target_data_pos['coexpression'])
coexpression_data_pos = coexpression_data_pos[~np.isnan(coexpression_data_pos)]
coexpression_data_neg = np.array(p_db.target_data_neg['coexpression'])
coexpression_data_neg = coexpression_data_neg[~np.isnan(coexpression_data_neg)]


pc = np.zeros(5)
percentiles = [16.666, 33.333, 50, 66.666, 83.333]
for idx, percentile in enumerate(percentiles):
    pc[idx] = np.percentile(coexpression_data,percentile)
    plt.axvline(pc[idx], color='k', linestyle='--')
plt.savefig('../course-project/coexpression.png', format='png', dpi=600)