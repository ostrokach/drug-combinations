import pandas as pd
import cPickle as pickle

import numpy as np
import matplotlib.pyplot as plt

def make_plot(mutation_type, mutation_stats, title_string):
    bar_labels = (" Can't map\n to uniprot ",
                  " No structural\n templates ", 
    #              " Interface\n templates ", 
                  " Only core\n templates ", 
                  " Interface\n and core\n templates ")
    n_groups = len(bar_labels)
    
    # Plot driver mutation statistics                       
    fig, ax = plt.subplots()
    
    index = np.arange(n_groups) + 0.35
    bar_width = 0.7
    opacity = 0.4
    error_config = {'ecolor': '0.3'}
    
    rects1 = plt.bar(index, mutation_stats, bar_width, alpha=opacity, color='b')
    

    plt.xticks(index + bar_width/2, bar_labels, size='large', rotation=0, ha='center')
    plt.xlim(0, 4.35)
    
    plt.ylabel('Number of %s' % binning_by, size='large')
    plt.title(title_string, size='large')    
    plt.tight_layout()
    plt.show()
    plt.savefig('/home/alexey/documents/presentations/subgroup-meetings/splicing/131204/%s-by-%s.png' % (mutation_type, binning_by), format='png', dpi=600)
    plt.close()

if False:
    driver_stats = categorise_mutations(drivers)
    make_plot('Driver', driver_stats, 'Driver mutations in COSMIC (> 5 records, > 3 pubmeds)')
    
    passenger_stats = categorise_mutations(passengers)
    make_plot('Passenger', passenger_stats, 'Passenger mutations in COSMIC (< 5 records, < 3 pubmeds)')



