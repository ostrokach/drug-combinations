# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 18:12:49 2013

@author: alexey
"""

def split_pos_neg(filename):
    extension= filename[-4:]
    
    ih = open(filename, 'r')
    oh_pos = open(filename[:-4] + '-pos' + extension, 'w')
    oh_neg = open(filename[:-4] + '-neg' + extension, 'w')
    
    for line_number, line in enumerate(ih):
        try:
            row = line.strip().split('\t')
            if row[0][:3] == 'Pos':
                oh_pos.writelines('\t'.join([row[1], row[2]]) + '\n')
            elif row[0][:3] == 'Neg':
                oh_neg.writelines('\t'.join([row[1], row[2]]) + '\n')
            else:
                print "Did not write this line:", line
        except IndexError:
            print filename
            print line_number, line
    ih.close
    oh_pos.close()
    oh_neg.close()


if __name__ == '__main__':
    
    filenames = [
        '/home/alexey/working/chemical-interactions/data/targets/Coexpression/Coexpression.txt',
        '/home/alexey/working/chemical-interactions/data/targets/Essentiality/Essentiality.txt',
        '/home/alexey/working/chemical-interactions/data/targets/FxnAssociation/FxnAssociation.txt',
        '/home/alexey/working/chemical-interactions/data/targets/GeneticInteraction/GeneticInteraction.txt',
        '/home/alexey/working/chemical-interactions/data/targets/GOAll_SharedFraction/GOAll_SharedFraction.txt',
        '/home/alexey/working/chemical-interactions/data/targets/GOBP_Pattern/GO_BP_Pattern.txt',
        '/home/alexey/working/chemical-interactions/data/targets/GOBP_SharedFraction/GOBP_SharedFraction.txt',
        '/home/alexey/working/chemical-interactions/data/targets/GOCC_Pattern/GO_CC_Pattern.txt',
        '/home/alexey/working/chemical-interactions/data/targets/GOCC_SharedFraction/GOCC_SharedFraction.txt',
        '/home/alexey/working/chemical-interactions/data/targets/GOMF_Pattern/GO_MF_Pattern.txt',
        '/home/alexey/working/chemical-interactions/data/targets/GOMF_SharedFraction/GOMF_SharedFraction.txt',
        '/home/alexey/working/chemical-interactions/data/targets/Pathway_Pattern/Pathway_Pattern.txt',
        '/home/alexey/working/chemical-interactions/data/targets/Pathway_SharedFraction/Pathway_SharedFraction.txt',
        '/home/alexey/working/chemical-interactions/data/targets/RegNetwork/RegNetwork.txt',
        '/home/alexey/working/chemical-interactions/data/targets/ShortestPath/ShortestPath.txt',
        ]
    
    filenames = [
    '/home/alexey/working/chemical-interactions/data/targets/ShortestPath/ShortestPathLength/CRG_ShortestPathLength.txt']
    
    for filename in filenames:
        split_pos_neg(filename)