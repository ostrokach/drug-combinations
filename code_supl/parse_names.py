# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 19:42:37 2014

@author: alexey
"""
import numpy as np
import pandas as pd
import argparse

from collections import OrderedDict

from sqlalchemy import create_engine
from sqlalchemy import Column, Index
from sqlalchemy import Integer, Float, String, Boolean, Text

from sqlalchemy.orm import sessionmaker, relationship, backref, scoped_session

from sqlalchemy.ext.declarative import declarative_base


###############################################################################

parser= argparse.ArgumentParser('target file index')
parser.add_argument('key_index', type=int)
args = parser.parse_args()

###############################################################################

from bisect import bisect_left

def binary_search(a, x, lo=0, hi=None):   # can't use a to specify default for hi
    hi = hi if hi is not None else len(a) # hi defaults to len(a)
    pos = bisect_left(a,x,lo,hi)          # find insertion position
    return (pos if pos != hi and a[pos] == x else -1) # don't walk off the end

def symmetrize(a):
    return a + a.T - np.diag(a.diagonal())

###############################################################################


# Global path
path_to_target_db = '/home/kimlab1/strokach/working/chemical-interactions/data/targets/'

# Local paths
paths = OrderedDict()
paths['coexpression'] = ('coexpression/GeneExp_AllHuman.txt.clip', 1,)
paths['essentiality'] = ('essentiality/Essentiality_AllHuman.txt', 1,)
paths['functional_association'] = ('functional_association/FxnAss_AllHuman.txt', 0,)
paths['genetic_interactions'] = ('genetic_interactions/GeneticInteraction_AllHuman.txt', 0,)
paths['go_all_shared_fraction'] = ('go_all_shared_fraction/Fraction_All.txt', 1,)
paths['go_bp_shared_fraction'] = ('go_bp_shared_fraction/Fraction_BP.txt', 1,)
paths['go_cc_shared_fraction'] = ('go_cc_shared_fraction/Fraction_CC.txt', 1,)
paths['go_mf_shared_fraction'] = ('go_mf_shared_fraction/Fraction_MF.txt', 1,)
paths['go_slim_shared_fraction'] = ('go_slim_shared_fraction/Fraction_Slim.txt', 1,)
paths['pathway_shared_fraction'] = ('pathway_shared_fraction/Pathway_Fraction_All.txt', 1,)
paths['shortest_path_length'] = ('shortest_path_length/ShortestPathLength_AllHuman.txt', 1,)


# Local paths
paths = OrderedDict()
paths['CRG_Betwee_AllHuman'] = ('CRG_Betwee_AllHuman.txt', 1, 2)
paths['CRG_Degree_AllHuman'] = ('CRG_Degree_AllHuman.txt', 1, 2)
paths['Fraction_All'] = ('Fraction_All.txt', 1, 2)
paths['Fractio_MF'] = ('Fraction_MF.txt', 1, 2)
paths['GeneExp_AllHuman'] = ('GeneExp_AllHuman.txt', 1, 2)
paths['Pathway_Fraction_All'] = ('Pathway_Fraction_All.txt', 1, 2)
paths['ShortestPathLength_AllHuman'] = ('ShortestPathLength_AllHuman.txt', 1, 2)
paths['CRG_Closen_AllHuman'] = ('CRG_Closen_AllHuman.txt', 1, 2)
paths['CRG_Neighbor_AllHuman'] = ('CRG_Neighbor_AllHuman.txt', 1, 2)
paths['Fraction_BP'] = ('Fraction_BP.txt', 1, 2)
paths['Fraction_Slim'] = ('Fraction_Slim.txt', 1, 2)
paths['GeneticInteraction_AllHuman'] = ('GeneticInteraction_AllHuman.txt', 0, None)
paths['PhyloPro_AllHuman'] = ('PhyloPro_AllHuman.txt', 1, 2)
paths['CRG_CluCoe_AllHuman'] = ('CRG_CluCoe_AllHuman.txt', 1, 2)
paths['Essentiality_AllHuman'] = ('Essentiality_AllHuman.txt', 1, 2)
paths['Fraction_CC'] = ('Fraction_CC.txt', 1, 2)
paths['FxnAss_AllHuman'] = ('FxnAss_AllHuman.txt', 0, 2)
paths['GetInt_ShortstPahtLength_AllHuman'] = ('GetInt_ShortstPahtLength_AllHuman.txt', 1, 2)
paths['RH_GetInt_FDR_0'] = ('RH_GetInt_FDR_0.txt', 0, 1)
paths['RegNet_AllHuman'] = ('RegNet_AllHuman.txt', 1, 2)


###############################################################################

Base = declarative_base()

sql_flavor = 'mysql'


if sql_flavor.split('_')[0] == 'sqlite':
    binary_collation = 'RTRIM' # same as binary, except that trailing space characters are ignored.
    string_collation = 'NOCASE'
elif sql_flavor.split('_')[0] == 'mysql':
    binary_collation = 'utf8_bin'
    string_collation = 'utf8_unicode_ci'
else:
    raise Exception('Unknown database type!')


class Targets(Base):
    __tablename__ = 'targets_2'

    ensp_1 = Column(Integer, primary_key=True)
    ensp_2 = Column(Integer, primary_key=True)

    CRG_Betwee_AllHuman = Column(Float)
    CRG_Closen_AllHuman = Column(Float)
    CRG_CluCoe_AllHuman = Column(Float)
    CRG_Degree_AllHuman = Column(Float)
    CRG_Neighbor_AllHuman = Column(Float)
    Essentiality_AllHuman = Column(Float)
    Fraction_All = Column(Float)
    Fraction_BP = Column(Float)
    Fraction_CC = Column(Float)
    Fraction_MF = Column(Float)
    Fraction_Slim = Column(Float)
    FxnAss_AllHuman = Column(Integer)
    GeneExp_AllHuman = Column(Float)
    GeneticInteraction_AllHuman = Column(Integer)
    GetInt_ShortstPahtLength_AllHuman = Column(Integer)
    Pathway_Fraction_All = Column(Float)
    PhyloPro_AllHuman = Column(Float)
    RegNet_AllHuman = Column(Float)
    RH_GetInt_FDR_0 = Column(Float)
    ShortestPathLength_AllHuman = Column(Integer)

    Index('ensp', 'ensp_1', 'ensp_2', unique=True)


class PositiveSet(Base):
    __tablename__ = 'positive_set'
    ensp_1 = Column(Integer, primary_key=True)
    ensp_2 = Column(Integer, primary_key=True)
    drug_id_1 = Column(String(15), nullable=False, primary_key=True)
    drug_id_2 = Column(String(15), nullable=False, primary_key=True)
    __table_args__ = (Index('ensp_pair_pos', 'ensp_1', 'ensp_2'),
                      Index('drug_id_pair_pos', 'drug_id_1', 'drug_id_2'),)


class NegativeSet(Base):
    __tablename__ = 'negative_set'
    ensp_1 = Column(Integer, primary_key=True)
    ensp_2 = Column(Integer, primary_key=True)
    drug_id_1 = Column(String(15), nullable=False, primary_key=True)
    drug_id_2 = Column(String(15), nullable=False, primary_key=True)
    __table_args__ = (Index('ensp_pair_neg', 'ensp_1', 'ensp_2'),
                      Index('drug_id_pair_neg', 'drug_id_1', 'drug_id_2'),)

################################################################################
## Run for every dataset to get an exhaustive list of target proteins

#with open(path_to_target_db + paths[key][0]) as fh:
#    if key == 'functional_association':
#        fh.next()
#    if True:
#        # Read the file in one go
#        ensp_pairs = np.loadtxt(fh, delimiter='\t', usecols=(paths[key][1],), dtype='S31')
#    else:
#        # Read the file one line at a time, stopping after a particular line
#        ensp_pairs = []
#        for counter, line in enumerate(fh):
#            if counter > 50:
#                break
#            ensp_pairs.append(line.split('\t')[paths[key][1]])
#
#set_of_names = set()
#
#for pair in ensp_pairs:
#    try:
#        c1, c2 = pair.split('_')
#        assert c1[:4] == 'ENSP' and c2[:4] == 'ENSP'
#        c1 = int(c1[4:])
#        c2 = int(c2[4:])
#        set_of_names.add(c1)
#        set_of_names.add(c2)
#    except ValueError:
#        print pair
#        raise
#
#del ensp_pairs
#
#with open(path_to_target_db + paths[key][0] + '.unique_names.pickle', 'w') as fh:
#    pickle.dump(set_of_names, fh)


###############################################################################
# Connect to the database
#engine = create_engine('mysql://root:kim630@192.168.6.19:3306/chemical_interactions')
engine = create_engine('postgresql://strokach:@192.168.6.19:5432/chemical_interactions')
Base.metadata.drop_all(engine)
Base.metadata.create_all(engine)
Session = sessionmaker()
Session.configure(bind=engine)
session = Session()

#
################################################################################
drug_pairs_filename='../data/Drug_Pairs_PositiveSet.txt'
drugs_targets_pos_filename='../data/Target_Pairs_Rel2_Positive.txt'
drugs_targets_neg_filename='../data/Target_Pairs_Rel2_Negative.txt'


# DataFrame to convert unique ids Clare assigned to positive drug pairs, to STITCH drug accessions
with open(drug_pairs_filename, 'r') as ih:
    drug_pairs = pd.read_csv(ih, sep='\t', skiprows=1, names=['drug_pair_id', 'drug_pair'], dtype=str)
    drug_pairs.set_index('drug_pair_id', drop=False, inplace=True, verify_integrity=False)

# DataFrame to link all drug_pairs with all target_pairs in the positive set
drugs_targets_pos = set()
with open(drugs_targets_pos_filename, 'r') as ih:
    for line in ih:
        row = [ l.strip() for l in line.split('\t') ]
        target_pair = row[0]
        drug_pair_ids = row[2].replace('[','').replace(']','').replace("'", '').split(', ')
        for drug_pair_id in drug_pair_ids:
            drug_pair = drug_pairs.loc[drug_pair_id]['drug_pair']
            if (drug_pair, target_pair,) in drugs_targets_pos:
                print (drug_pair, target_pair,), 'already in positive set!'
                continue
            drugs_targets_pos.add((drug_pair, target_pair,))
drugs_targets_pos = pd.DataFrame(list(drugs_targets_pos),  columns=['drug_pair', 'target_pair'])

for line_number, (idx, df) in enumerate(drugs_targets_pos.iterrows()):
    ensps = df['target_pair'].split('_')
    assert (ensps[0][:4] == 'ENSP') & (ensps[1][:4] == 'ENSP')
    ensps = [int(e[4:]) for e in ensps]
    drugs = df['drug_pair'].split('_')
    if drugs[0][:3] == 'CID':
        drugs[0] = str(int(drugs[0][4:]))
    if drugs[1][:3] == 'CID':
        drugs[1] = str(int(drugs[1][4:]))
    if ensps[1] < ensps[0]:
        ensps[0], ensps[1] = ensps[1], ensps[0]
        drugs[0], drugs[1] = drugs[1], drugs[0]
    rec = PositiveSet()
    rec.ensp_1 = ensps[0]
    rec.ensp_2 = ensps[1]
    rec.drug_id_1 = drugs[0]
    rec.drug_id_2 = drugs[1]
    session.add(rec)
    if line_number % 1000 == 0:
        session.flush()
        print line_number
session.commit()


# DataFrame to link all drug_pairs with all target_pairs in the negative set
drugs_targets_neg = set()
with open(drugs_targets_neg_filename, 'r') as ih:
    for line in ih:
        row = [ l.strip() for l in line.split('\t') ]
        target_pair = row[0]
        drug_pairs = row[2].replace('[','').replace(']','').replace("'", '').split(', ')
        for drug_pair in drug_pairs:
            if (drug_pair, target_pair,) in drugs_targets_pos:
                print (drug_pair, target_pair,), 'are also in positive set!'
                continue
            if (drug_pair, target_pair,) in drugs_targets_neg:
                print (drug_pair, target_pair,), 'already in negative set!'
                continue
            drugs_targets_neg.add((drug_pair, target_pair,))
drugs_targets_neg = pd.DataFrame(list(drugs_targets_neg), columns=['drug_pair', 'target_pair'])

for line_number, (idx, df) in enumerate(drugs_targets_neg.iterrows()):
    ensps = df['target_pair'].split('_')
    assert (ensps[0][:4] == 'ENSP') & (ensps[1][:4] == 'ENSP')
    ensps = [int(e[4:]) for e in ensps]
    drugs = df['drug_pair'].split('_')
    if drugs[0][:3] == 'CID':
        drugs[0] = str(int(drugs[0][4:]))
    if drugs[1][:3] == 'CID':
        drugs[1] = str(int(drugs[1][4:]))
    if ensps[1] < ensps[0]:
        ensps[0], ensps[1] = ensps[1], ensps[0]
        drugs[0], drugs[1] = drugs[1], drugs[0]
    rec = NegativeSet()
    rec.ensp_1 = ensps[0]
    rec.ensp_2 = ensps[1]
    rec.drug_id_1 = drugs[0]
    rec.drug_id_2 = drugs[1]
    session.add(rec)
    if line_number % 1000 == 0:
        session.flush()
        print line_number
session.commit()
