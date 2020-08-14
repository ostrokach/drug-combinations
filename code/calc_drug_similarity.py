# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 00:00:26 2013

@author: alexey
"""

from rdkit import Chem, DataStructs

def get_similarity_all(fp1, fp2):
    """
    Get similarity score for fingerprints that are supplied always as SparseBitVect
    RDKit has the following similarity measures:
        Tanimoto, Dice, Cosine, Sokal, Russel, Kulczynski, McConnaughey, and Tversky.
    """
    similarity_scores = [
        DataStructs.TanimotoSimilarity(fp1,fp2),
        DataStructs.DiceSimilarity(fp1,fp2),
        DataStructs.CosineSimilarity(fp1,fp2),
#        DataStructs.SokalSimilarity(fp1,fp2),
        DataStructs.RusselSimilarity(fp1,fp2),
        DataStructs.KulczynskiSimilarity(fp1,fp2),
        DataStructs.McConnaugheySimilarity(fp1,fp2)]

    return similarity_scores


def get_similarity_subset(fp1, fp2):
    """
    Get similarity score for fingerprints that are supplied as ExplicitBitVect
    or some other format.
    The following similarity metrics work with different intput formats:
        Tanimoto, Dice
    """
    similarity_scores = [
        DataStructs.TanimotoSimilarity(fp1,fp2),
        DataStructs.DiceSimilarity(fp1,fp2)]

    return similarity_scores


def calculate_similarity_vector(smile_pair):
    """
    Calculate fingerprints between two smile terms using different fingerprinters,
    and use different similarity metrics to calculate the difference between those fingerprints.
    """
#    smile1, smile2 = smile_pair.split('_')
    smile1, smile2 = smile_pair

    mol1 = Chem.MolFromSmiles(smile1)
    mol2 = Chem.MolFromSmiles(smile2)

    molecule_similarity = list()

    # RDK topological fingerprint for a molecule
    fp1 = Chem.RDKFingerprint(mol1)
    fp2 = Chem.RDKFingerprint(mol2)
    molecule_similarity.extend(get_similarity_all(fp1, fp2))
    #print 'RDK fingerprint: ', DataStructs.KulczynskiSimilarity(fp1,fp2)

    ## LayeredFingerprint, a fingerprint using SMARTS patterns
    #fp1 = Chem.LayeredFingerprint(mol1)
    #fp2 = Chem.LayeredFingerprint(mol2)
    #print 'RDK fingerprint: ', DataStructs.TanimotoSimilarity(fp1,fp2)

    # PatternFingerprint, a fingerprint using SMARTS patterns
    #fp1 = Chem.PatternFingerprint(mol1)
    #fp2 = Chem.PatternFingerprint(mol2)
    #print 'RDK fingerprint: ', DataStructs.TanimotoSimilarity(fp1,fp2)

    ###############################################################################

    # Topological Fingerprints
    # Uses Chem.RDKFingerprint internally, but with different parameters, I guess...
    # http://www.rdkit.org/docs/GettingStartedInPython.html#topological-fingerprints
    from rdkit.Chem.Fingerprints import FingerprintMols
    fp1 = FingerprintMols.FingerprintMol(mol1)
    fp2 = FingerprintMols.FingerprintMol(mol2)
    molecule_similarity.extend(get_similarity_all(fp1, fp2))
    #print 'RDK fingerprint: ', DataStructs.TanimotoSimilarity(fp1,fp2)

    ###############################################################################

    # MACCS Keys
    # There is a SMARTS-based implementation of the 166 public MACCS keys.
    # http://www.rdkit.org/docs/GettingStartedInPython.html#maccs-keys
    from rdkit.Chem import MACCSkeys
    fp1 = MACCSkeys.GenMACCSKeys(mol1)
    fp2 = MACCSkeys.GenMACCSKeys(mol2)
    molecule_similarity.extend(get_similarity_all(fp1, fp2))
    #print "RDK fingerprint: ", DataStructs.TanimotoSimilarity(fp1,fp2)

    ###############################################################################

    # Atom Pairs and Topological Torsions
    # Atom-pair descriptors [3] are available in several different forms.
    # The standard form is as fingerprint including counts for each bit instead of just zeros and ones:
    # http://www.rdkit.org/docs/GettingStartedInPython.html#atom-pairs-and-topological-torsions
    from rdkit.Chem.AtomPairs import Pairs
    fp1 = Pairs.GetAtomPairFingerprintAsBitVect(mol1)
    fp2 = Pairs.GetAtomPairFingerprintAsBitVect(mol2)
    molecule_similarity.extend(get_similarity_all(fp1, fp2))
    #print "RDK fingerprint: ", DataStructs.DiceSimilarity(fp1,fp2)
    from rdkit.Chem.AtomPairs import Torsions
    fp1 = Torsions.GetTopologicalTorsionFingerprint(mol1)
    fp2 = Torsions.GetTopologicalTorsionFingerprint(mol2)
    molecule_similarity.extend(get_similarity_subset(fp1, fp2))
    #print "RDK fingerprint: ", DataStructs.TanimotoSimilarity(fp1,fp2)

    ###############################################################################

    # Morgan Fingerprints (Circular Fingerprints)
    #This family of fingerprints, better known as circular fingerprints [5],
    #is built by applying the Morgan algorithm to a set of user-supplied atom invariants.
    #When generating Morgan fingerprints, the radius of the fingerprint must also be provided...
    # http://www.rdkit.org/docs/GettingStartedInPython.html#morgan-fingerprints-circular-fingerprints
    from rdkit.Chem import rdMolDescriptors
    fp1 = rdMolDescriptors.GetMorganFingerprint(mol1,2)
    fp2 = rdMolDescriptors.GetMorganFingerprint(mol2,2)
    molecule_similarity.extend(get_similarity_subset(fp1, fp2))

    fp1 = rdMolDescriptors.GetMorganFingerprint(mol1,2,useFeatures=True)
    fp2 = rdMolDescriptors.GetMorganFingerprint(mol2,2,useFeatures=True)
    molecule_similarity.extend(get_similarity_subset(fp1, fp2))

    #print "RDK fingerprint: ", DataStructs.TanimotoSimilarity(fp1,fp2)

    ###############################################################################

    return molecule_similarity


if __name__ == '__main__':
    import pandas as pd
    import argparse

    if False:
        # Old code for parsing Clare's predictor data
        input_path = '/home/kimlab1/strokach/working/chemical_interactions/drug_similarity/'
        output_path = '/home/kimlab1/strokach/working/chemical_interactions/drug_similarity/'
        input_filename = 'All_TestedDrugs.txt'
        input_file = pd.read_csv(input_path + input_filename, sep='\t')
        input_file_subset = input_file[['STITCH_PubChem(CID)','SMILE']]
        del input_file
        input_file_subset['cid'] = input_file_subset['STITCH_PubChem(CID)'].apply(lambda x: int(x[3:]))
        input_file_subset['smile'] = input_file_subset['SMILE']
        del input_file_subset['STITCH_PubChem(CID)']
        del input_file_subset['SMILE']
        input_file_subset.sort(['cid'], inplace=True)

        predictor_filenames = [
            '/home/kimlab1/cjeon/DrugCombination_2/IndependentValidationSet_PredictorI.txt',
            '/home/kimlab1/cjeon/DrugCombination_2/PredictorI.txt',
            '/home/kimlab1/cjeon/DrugCombination_2/PredictorI_HighScored.txt']

        def split_and_sort(drug_pair):
            a, b = [int(x[3:]) for x in drug_pair.split('_')]
            if a <= b:
                return [a, b]
            else:
                return [b, a]

        predictor_1 = pd.read_csv(predictor_filenames[1], sep='\t')
        predictor_1['cid_a'] = predictor_1['DrugPairs'].apply(lambda x: split_and_sort(x)[0])
        predictor_1['cid_b'] = predictor_1['DrugPairs'].apply(lambda x: split_and_sort(x)[1])

        predictor_files = []
        for predictor_filename in predictor_filenames:
            predictor_files.append(pd.read_csv(predictor_filename, sep='\t'))
        predictor_file = pd.concat(predictor_files)
        predictor_file['cid_a'] = predictor_file['DrugPairs'].apply(lambda x: split_and_sort(x)[0])
        predictor_file['cid_b'] = predictor_file['DrugPairs'].apply(lambda x: split_and_sort(x)[0])

        predictor_file_subset = predictor_file[['cid_a', 'cid_b', 'DrugPairs']]
#        del predictor_file
        predictor_file_subset['cid_a'] = predictor_file_subset['DrugPairs'].apply(lambda x: split_and_sort(x)[0])
        predictor_file_subset['cid_b'] = predictor_file_subset['DrugPairs'].apply(lambda x: split_and_sort(x)[1])
        predictor_file_subset.drop_duplicates(inplace=True)
        predictor_file_subset = predictor_file_subset \
            .merge(input_file_subset, how='left', left_on='cid_a', right_on='cid') \
            .merge(input_file_subset, how='left', left_on='cid_b', right_on='cid', suffixes=('_1', '_2'))
        assert sum([cid_a!=cid_1 for (cid_a, cid_1) in predictor_file_subset[['cid_a', 'cid_1']].values]) == 0
        assert sum([cid_b!=cid_2 for (cid_b, cid_2) in predictor_file_subset[['cid_b', 'cid_2']].values]) == 0
#        del predictor_file_subset['cid_a']
#        del predictor_file_subset['cid_b']
    elif False:
        # New code for getting chemical-chemical similarity for the entire database
        parser = argparse.ArgumentParser()
        parser.add_argument('--bin_idx', type=int, default=0)
        args = parser.parse_args()

        path_to_data = '/home/kimlab1/strokach/working/databases/chemical_interactions/version_2.1/'
        input_filename = 'cid_smile_pairs_remaining.tsv'
        output_filename = input_filename + '.w_chem_similarity.idx_{}'.format(args.bin_idx)

        cid_pair_df = pd.read_csv(path_to_data + input_filename, sep='\t')

        no_bins = 1
        bin_size = len(cid_pair_df) / no_bins + 1
        from_idx = bin_size * args.bin_idx
        to_idx = min([bin_size * (args.bin_idx + 1), len(cid_pair_df)])

        similarity_metric_all = ['Tanimoto', 'Dice', 'Cosine', 'Russel', 'Kulczynski', 'McConnaughey']
        similarity_metric_subset = ['Tanimoto', 'Dice']
        column_names = (
            ['RDKFingerprint_' + sim for sim in similarity_metric_all] +
            ['FingerprintMol_' + sim for sim in similarity_metric_all] +
            ['MACCSkeys_' + sim for sim in similarity_metric_all] +
            ['AtomPairFingerprint_' + sim for sim in similarity_metric_all] +
            ['TopologicalTorsionFingerprint_' + sim for sim in similarity_metric_subset] +
            ['MorganFingerprintR2_' + sim for sim in similarity_metric_subset] +
            ['MorganFingerprintR2withFeatures_' + sim for sim in similarity_metric_subset])

        cid_pair_df_subset = cid_pair_df[from_idx:to_idx].copy()
        chem_sim_matrix = [calculate_similarity_vector([x[0],x[1]]) for x in
            cid_pair_df_subset[['smile_1','smile_2']].values]
        chem_sim_df = pd.DataFrame(chem_sim_matrix, columns=column_names)
        chem_sim_df.index = cid_pair_df_subset.index
        cid_pair_chem_sim_df = cid_pair_df_subset.join(chem_sim_df)
        del cid_pair_chem_sim_df['smile_1']
        del cid_pair_chem_sim_df['smile_2']
        cid_pair_chem_sim_df.to_csv(path_to_data + output_filename, sep='\t', index=False)









#if __name__ == '__main__':
#    import pandas as pd
#    import psycopg2
#    import cPickle as pickle
#
#    predictor_1_filename = 'PredictorI.txt'
#    predictor_2_filename = 'PredictorII.txt'
#    conn = psycopg2.connect(database='kimlab', user='postgres', host='192.168.6.19', port='5432')
#
#    set_of_bad_smiles = set(['[HH]', '[O-]Cl(=O)=O.[Na+]'])
#
#    def get_smile_pair(drug_pair):
#        def get_smile(cid):
#            curs = conn.cursor()
#            command = "SELECT smiles_string FROM stitch.chemicals WHERE chemical = '%s'" % cid
#            curs.execute(command)
#            smiles = curs.fetchall()
#            if smiles:
#                return smiles[0][0]
#            else:
#                return ''
#        cid_1, cid_2 = drug_pair.split('_')
#        smile_1 = get_smile(cid_1)
#        smile_2 = get_smile(cid_2)
#        return smile_1 + '_' + smile_2
#
#    path_to_data = '/home/kimlab1/strokach/working/chemical_interactions/data/version_2/parsed_outputs/'
#    #similarity_metric_all = ['Tanimoto', 'Dice', 'Cosine', 'Sokal', 'Russel', 'Kulczynski', 'McConnaughey']
#    similarity_metric_all = ['Tanimoto', 'Dice', 'Cosine', 'Russel', 'Kulczynski', 'McConnaughey']
#    similarity_metric_subset = ['Tanimoto', 'Dice']
#    column_names = (
#        ['RDKFingerprint_' + sim for sim in similarity_metric_all] +
#        ['FingerprintMol_' + sim for sim in similarity_metric_all] +
#        ['MACCSkeys_' + sim for sim in similarity_metric_all] +
#        ['AtomPairFingerprint_' + sim for sim in similarity_metric_all] +
#        ['TopologicalTorsionFingerprint_' + sim for sim in similarity_metric_subset] +
#        ['MorganFingerprintR2_' + sim for sim in similarity_metric_subset] +
#        ['MorganFingerprintR2withFeatures_' + sim for sim in similarity_metric_subset])
#
#    if False:
#        pred_1_df = pd.read_csv(path_to_data + predictor_1_filename, sep='\t')
#        pred_1_df['smile_codes'] = pred_1_df['DrugPair'].apply(get_smile_pair)
#        pred_1_df.to_csv(path_to_data + predictor_1_filename.replace('.txt','.with_smiles.txt'), index=False, sep='\t')
#    else:
#        pred_1_df = pd.read_csv(path_to_data + predictor_1_filename.replace('.txt','.with_smiles.txt'), sep='\t')
#    pred_1_df['Type'] = pred_1_df['Type'].apply(lambda x: x.strip())
#    pred_1_df['DrugPair'] = pred_1_df['DrugPair'].apply(lambda x: x.strip())
#    pred_1_df['TargetPair'] = pred_1_df['TargetPair'].apply(lambda x: x.strip())
#
#    if True:
#        pred_1_features = pred_1_df['smile_codes'].apply(calculate_similarity_vector)
#        with open(path_to_data + predictor_1_filename.replace('.txt','.chemical_similarity.pickle'), 'wb') as fh:
#            pickle.dump(pred_1_features, fh)
#    else:
#        with open(path_to_data + predictor_1_filename.replace('.txt','.chemical_similarity.pickle'), 'rb') as fh:
#            pred_1_features = pickle.load(fh)
#    pred_1_features_df = pd.DataFrame(pred_1_features.values.tolist(), columns=column_names)
#    pred_1_df_all = pred_1_df.join(pred_1_features_df)
#    del pred_1_df_all['smile_codes']
#    pred_1_atc_sim = pd.read_csv(path_to_data + 'drugs/ATC_PredictorI.txt', sep='\t')
#    pred_1_atc_sim['Type'] = pred_1_atc_sim['Type'].apply(lambda x: x.strip())
#    pred_1_atc_sim['DrugPair'] = pred_1_atc_sim['DrugPair'].apply(lambda x: x.strip())
#    pred_1_df_all = pred_1_df_all.merge(pred_1_atc_sim, how='left', on=['Type', 'DrugPair'])
#
#    pred_1_side_effect_sim = pd.read_csv(path_to_data + 'drugs/SideEffect_PredictorI.txt', sep='\t')
#    pred_1_side_effect_sim['Type'] = pred_1_side_effect_sim['Type'].apply(lambda x: x.strip())
#    pred_1_side_effect_sim['DrugPair'] = pred_1_side_effect_sim['DrugPair'].apply(lambda x: x.strip())
#    pred_1_df_all = pred_1_df_all.merge(pred_1_side_effect_sim, how='left', on=['Type', 'DrugPair'])
#
#    pred_1_df_all.to_csv(path_to_data + predictor_1_filename.replace('.txt','.with_drug_features.txt'), index=False, sep='\t')
#
#
#    if False:
#        pred_2_df = pd.read_csv(path_to_data + predictor_2_filename, sep='\t')
#        pred_2_df['smile_codes'] = pred_2_df['DrugPair'].apply(get_smile_pair)
#        pred_2_df.to_csv(path_to_data + predictor_2_filename.replace('.txt','.with_smiles.txt'), index=False, sep='\t')
#    else:
#        pred_2_df = pd.read_csv(path_to_data + predictor_2_filename.replace('.txt','.with_smiles.txt'), sep='\t')
#    pred_2_df['Type'] = pred_2_df['Type'].apply(lambda x: x.strip())
#    pred_2_df['DrugPair'] = pred_2_df['DrugPair'].apply(lambda x: x.strip())
#    pred_2_df['TargetPair'] = pred_2_df['TargetPair'].apply(lambda x: x.strip())
#
#    if True:
#        pred_2_features = pred_2_df['smile_codes'].apply(calculate_similarity_vector)
#        with open(path_to_data + predictor_2_filename.replace('.txt','.chemical_similarity.pickle'), 'wb') as fh:
#            pickle.dump(pred_2_features, fh)
#    else:
#        with open(path_to_data + predictor_2_filename.replace('.txt','.chemical_similarity.pickle'), 'rb') as fh:
#            pred_2_features = pickle.load(fh)
#    pred_2_features_df = pd.DataFrame(pred_2_features.values.tolist(), columns=column_names)
#    pred_2_df_all = pred_2_df.join(pred_2_features_df)
#    del pred_2_df_all['smile_codes']
#    pred_2_df_all.to_csv(path_to_data + predictor_2_filename.replace('.txt','.with_drug_features.txt'), index=False, sep='\t')
#
#
#    ###############################################################################
#    # Legacy code
#    if False:
#        if (os.path.isfile(path_to_data + 'pickled/chemical_similarity_pos' + '.pickle')
#        and os.path.isfile(path_to_data + 'pickled/chemical_similarity_neg' + '.pickle')):
#            # Load pickled data from a previous session if it exists
#            with open(path_to_data + 'pickled/smiles_terms_pos' + '.pickle', 'rb') as fh:
#                smiles_terms_pos = pickle.load(fh)
#            with open(path_to_data + 'pickled/chemical_similarity_pos' + '.pickle', 'rb') as fh:
#                chemical_similarity_pos = pickle.load(fh)
#            with open(path_to_data + 'pickled/smiles_terms_neg' + '.pickle', 'rb') as fh:
#                smiles_terms_neg = pickle.load(fh)
#            with open(path_to_data + 'pickled/chemical_similarity_neg' + '.pickle', 'rb') as fh:
#                chemical_similarity_neg = pickle.load(fh)
#        else:
#            # Calculate data from scratch
#            drug_pairs_index_pos = [ l.split('_') for l in project_db.drug_data_pos.index ]
#            smiles_terms_pos, chemical_similarity_pos = fill_sim_data(drug_pairs_index_pos)
#
#            drug_pairs_index_neg = [ l.split('_') for l in project_db.drug_data_neg.index ]
#            smiles_terms_neg, chemical_similarity_neg = fill_sim_data(drug_pairs_index_neg)
#
#            # Save calculated data
#            with open(path_to_data + 'pickled/smiles_terms_pos' + '.pickle', 'wb') as fh:
#                pickle.dump(smiles_terms_pos, fh)
#            with open(path_to_data + 'pickled/chemical_similarity_pos' + '.pickle', 'wb') as fh:
#                pickle.dump(chemical_similarity_pos, fh)
#            with open(path_to_data + 'pickled/smiles_terms_neg' + '.pickle', 'wb') as fh:
#                pickle.dump(smiles_terms_neg, fh)
#            with open(path_to_data + 'pickled/chemical_similarity_neg' + '.pickle', 'wb') as fh:
#                pickle.dump(chemical_similarity_neg, fh)
#
#
#        similarity_metric_all = ['Tanimoto', 'Dice', 'Cosine', 'Sokal', 'Russel', 'Kulczynski', 'McConnaughey']
#        similarity_metric_subset = ['Tanimoto', 'Dice']
#        column_names = (
#            ['RDKFingerprint-' + sim for sim in similarity_metric_all] +
#            ['FingerprintMol-' + sim for sim in similarity_metric_all] +
#            ['MACCSkeys-' + sim for sim in similarity_metric_all] +
#            ['AtomPairFingerprint-' + sim for sim in similarity_metric_all] +
#            ['TopologicalTorsionFingerprint-' + sim for sim in similarity_metric_subset] +
#            ['MorganFingerprintR2-' + sim for sim in similarity_metric_subset] +
#            ['MorganFingerprintR2withFeatures-' + sim for sim in similarity_metric_subset])
#
#
#        chemical_similarity_pos_base_df = pd.DataFrame(data=smiles_terms_pos,
#                                                       index=project_db.drug_data_pos.index,
#                                                       columns=['SMILES_drug1', 'SMILES_drug2'])
#        chemical_similarity_neg_base_df = pd.DataFrame(data=smiles_terms_neg,
#                                                       index=project_db.drug_data_neg.index,
#                                                       columns=['SMILES_drug1', 'SMILES_drug2'])
#        chemical_similarity_pos_base_df.to_csv(path_to_data + 'drugs/chemical-similarity/chemical_similarity_base_pos.tsv', sep='\t')
#        chemical_similarity_neg_base_df.to_csv(path_to_data + 'drugs/chemical-similarity/chemical_similarity_base_neg.tsv', sep='\t')
#
#
#        chemical_similarity_pos_data_df = pd.DataFrame(data=chemical_similarity_pos,
#                                                  index=project_db.drug_data_pos.index,
#                                                  columns=column_names)
#        chemical_similarity_neg_data_df = pd.DataFrame(data=chemical_similarity_neg,
#                                                  index=project_db.drug_data_neg.index,
#                                                  columns=column_names)
#        chemical_similarity_pos_data_df.to_csv(path_to_data + 'drugs/chemical-similarity/chemical_similarity_data_pos.tsv', sep='\t')
#        chemical_similarity_neg_data_df.to_csv(path_to_data + 'drugs/chemical-similarity/chemical_similarity_data_neg.tsv', sep='\t')
#
#
#        chemical_similarity_pos_df = pd.DataFrame(data=np.concatenate((smiles_terms_pos, chemical_similarity_pos),axis=1),
#                                                  index=project_db.drug_data_pos.index,
#                                                  columns=['SMILES_drug1', 'SMILES_drug2']+column_names)
#        chemical_similarity_neg_df = pd.DataFrame(data=np.concatenate((smiles_terms_neg, chemical_similarity_neg),axis=1),
#                                                  index=project_db.drug_data_neg.index,
#                                                  columns=['SMILES_drug1', 'SMILES_drug2']+column_names)
#        chemical_similarity_pos_df.to_csv(path_to_data + 'drugs/chemical-similarity/chemical_similarity_all_pos.tsv', sep='\t')
#        chemical_similarity_neg_df.to_csv(path_to_data + 'drugs/chemical-similarity/chemical_similarity_all_neg.tsv', sep='\t')
#
#
