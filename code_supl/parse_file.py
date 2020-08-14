filenames = [
    '/tmp/DrugCombination_2/BioGridTopo_EB/BioGridTopo_EB_PredictorI.txt',
    '/tmp/DrugCombination_2/BioGridTopo_EB/BioGridTopo_EB_PredictorII.txt',
    '/tmp/DrugCombination_2/GetInt_EB/GetInt_EB_PredictorI.txt',
    '/tmp/DrugCombination_2/GetInt_EB/GetInt_EB_PredictorII.txt',
    '/tmp/DrugCombination_2/STRINGTopo_EB/STRINGTopo_EB_PredictorI.txt',
    '/tmp/DrugCombination_2/STRINGTopo_EB/STRINGTopo_EB_PredictorII.txt']
for filename in filenames:
    with open(filename) as ifh:
        ifh.readline()
        with open(filename.replace('.txt', '_parsed.txt'), 'w') as ofh:
            for line in ifh:
                row = line.strip().split('\t')
                a = row[0].strip()
                b = int(row[1].strip().split('_')[0][4:])
                c = int(row[1].strip().split('_')[1][4:])
                d = row[2].strip()
                e = row[3].strip()
                f = row[4].strip()
                g = row[5].strip()
                if b <= c:
                    ofh.writelines('\t'.join([a, str(b), str(c), d, e, f, g]) + '\n')
                else:
                    ofh.writelines('\t'.join([a, str(c), str(b), d, e, f, g]) + '\n')













filenames = [
'/tmp/DrugCombination_2/PredictorI.txt',
'/tmp/DrugCombination_2/PredictorI_HighScored.txt',
'/tmp/DrugCombination_2/PredictorII.txt',
'/tmp/DrugCombination_2/PredictorII_HighScored.txt']

for filename in filenames:
    with open(filename) as ifh:
        ifh.readline()
        with open(filename.replace('.txt', '_parsed.txt'), 'w') as ofh:
            for line in ifh:
                row = line.strip().split('\t')
                a = row[0].strip()
                b = int(row[1].strip().split('_')[0][4:])
                c = int(row[1].strip().split('_')[1][4:])
                d = int(row[2].strip().split('_')[0][3:])
                e = int(row[2].strip().split('_')[1][3:])
                if b <= c and d <= e:
                    ofh.writelines('\t'.join([a, str(b), str(c), str(d), str(e)]) + '\n')
                elif b <= c and d > e:
                    ofh.writelines('\t'.join([a, str(b), str(c), str(e), str(d)]) + '\n')
                elif b > c and d <= e:
                    ofh.writelines('\t'.join([a, str(c), str(b), str(d), str(e)]) + '\n')
                else:
                    ofh.writelines('\t'.join([a, str(c), str(b), str(e), str(d)]) + '\n')




### Phylo
filename = '/tmp/DrugCombination_2/Phylo/Phylo_AllHuman.txt'
with open(filename) as ifh:
    with open(filename.replace('.txt', '_parsed.txt'), 'w') as ofh:
        for line in ifh:
            row = line.strip().split('\t')
            a = int(row[0].strip().split('_')[0][4:])
            b = int(row[0].strip().split('_')[1][4:])
            c = row[1].strip()
            if a <= b:
                ofh.writelines('\t'.join([str(a), str(b), c]) + '\n')
            else:
                ofh.writelines('\t'.join([str(b), str(a), c]) + '\n')


### GO
filenames = [
    '/tmp/DrugCombination_2/GO/GO_AllPro_All.txt',
    '/tmp/DrugCombination_2/GO/GO_AllPro_BP.txt',
    '/tmp/DrugCombination_2/GO/GO_AllPro_CC.txt',
    '/tmp/DrugCombination_2/GO/GO_AllPro_MF.txt',
    '/tmp/DrugCombination_2/GO/GO_AllPro_Slim.txt',]

for filename in filenames:
    with open(filename) as ifh:
        with open(filename.replace('.txt', '_parsed.txt'), 'w') as ofh:
            for line in ifh:
                row = line.strip().split('\t')
                a = int(row[0].strip().split('_')[0][4:])
                b = int(row[0].strip().split('_')[1][4:])
                c = row[1].strip()
                if a <= b:
                    ofh.writelines('\t'.join([str(a), str(b), c]) + '\n')
                else:
                    ofh.writelines('\t'.join([str(b), str(a), c]) + '\n')


#~ ########################################################################


### GetInt

filename = '/tmp/DrugCombination_2/GetInt/GetIntTopo_Allhuman.txt'
with open(filename) as ifh:
    with open(filename.replace('.txt', '_parsed.txt'), 'w') as ofh:
        for line in ifh:
            row = line.strip().split('\t')
            a = row[0].strip()
            b = int(row[1].strip().split('_')[0][4:])
            c = int(row[1].strip().split('_')[1][4:])
            d = row[2].strip()
            e = row[3].strip()
            f = row[4].strip()
            g = row[5].strip()
            h = row[6].strip()
            i = row[7].strip()
            if b <= c:
                ofh.writelines('\t'.join([a, str(b), str(c), d, e, f, g, h, i]) + '\n')
            else:
                ofh.writelines('\t'.join([a, str(c), str(b), d, e, f, g, h, i]) + '\n')


filename = '/tmp/DrugCombination_2/GetInt/GetInt_EB_Allhuman.txt'
with open(filename) as ifh:
    with open(filename.replace('.txt', '_parsed.txt'), 'w') as ofh:
        for line in ifh:
            row = line.strip().split('\t')
            a = row[0].strip()
            b = int(row[1].strip().split('_')[0][4:])
            c = int(row[1].strip().split('_')[1][4:])
            d = row[2].strip()
            e = row[3].strip()
            f = row[4].strip()
            g = row[5].strip()
            if b <= c:
                ofh.writelines('\t'.join([a, str(b), str(c), d, e, f, g]) + '\n')
            else:
                ofh.writelines('\t'.join([a, str(c), str(b), d, e, f, g]) + '\n')


filename = '/tmp/DrugCombination_2/GetInt/GetIntTopo_NSP_Allhuman.txt'
with open(filename) as ifh:
    with open(filename.replace('.txt', '_parsed.txt'), 'w') as ofh:
        for line in ifh:
            row = line.strip().split('\t')
            a = row[0].strip()
            b = int(row[1].strip().split('_')[0][4:])
            c = int(row[1].strip().split('_')[1][4:])
            d = row[2].strip()
            if b <= c:
                ofh.writelines('\t'.join([a, str(b), str(c), d]) + '\n')
            else:
                ofh.writelines('\t'.join([a, str(c), str(b), d]) + '\n')



### PPI
filename = '/tmp/DrugCombination_2/PPI/BioGridTopo_Allhuman.txt'
with open(filename) as ifh:
    with open(filename.replace('.txt', '_parsed.txt'), 'w') as ofh:
        for line in ifh:
            row = line.strip().split('\t')
            a = row[0].strip()
            b = int(row[1].strip().split('_')[0][4:])
            c = int(row[1].strip().split('_')[1][4:])
            d = row[2].strip()
            e = row[3].strip()
            f = row[4].strip()
            g = row[5].strip()
            h = row[6].strip()
            i = row[7].strip()
            if b <= c:
                ofh.writelines('\t'.join([a, str(b), str(c), d, e, f, g, h, i]) + '\n')
            else:
                ofh.writelines('\t'.join([a, str(c), str(b), d, e, f, g, h, i]) + '\n')


filename = '/tmp/DrugCombination_2/PPI/BioGridTopo_EB_Allhuman.txt'
with open(filename) as ifh:
    with open(filename.replace('.txt', '_parsed.txt'), 'w') as ofh:
        for line in ifh:
            row = line.strip().split('\t')
            try:
                a = row[0].strip()
                b = int(row[1].strip().split('_')[0][4:])
                c = int(row[1].strip().split('_')[1][4:])
                d = row[2].strip()
                e = row[3].strip()
                f = row[4].strip()
                g = row[5].strip()
            except (IndexError, ValueError):
                print line
                print row
                continue
            if b <= c:
                ofh.writelines('\t'.join([a, str(b), str(c), d, e, f, g]) + '\n')
            else:
                ofh.writelines('\t'.join([a, str(c), str(b), d, e, f, g]) + '\n')


filename = '/tmp/DrugCombination_2/PPI/BioGridTopo_NSP_Allhuman.txt'
with open(filename) as ifh:
    with open(filename.replace('.txt', '_parsed.txt'), 'w') as ofh:
        for line in ifh:
            row = line.strip().split('\t')
            a = row[0].strip()
            b = int(row[1].strip().split('_')[0][4:])
            c = int(row[1].strip().split('_')[1][4:])
            d = row[2].strip()
            if b <= c:
                ofh.writelines('\t'.join([a, str(b), str(c), d]) + '\n')
            else:
                ofh.writelines('\t'.join([a, str(c), str(b), d]) + '\n')



########################################################################


#~ ### GeneExp
filename = '/tmp/DrugCombination_2/GeneExp/GeneExp_AllHuman.txt'
with open(filename) as ifh:
    with open(filename.replace('.txt', '_parsed.txt'), 'w') as ofh:
        for line in ifh:
            row = line.strip().split('\t')
            a = int(row[0].strip().split('_')[0][4:])
            b = int(row[0].strip().split('_')[1][4:])
            c = row[1].strip()
            if a <= b:
                ofh.writelines('\t'.join([str(a), str(b), c]) + '\n')
            else:
                ofh.writelines('\t'.join([str(b), str(a), c]) + '\n')


#~ ### GeneEss
filename = '/tmp/DrugCombination_2/GeneEss/Human_GeneEss_ENSP.txt'
with open(filename) as ifh:
    with open(filename.replace('.txt', '_parsed.txt'), 'w') as ofh:
        for line in ifh:
            row = line.strip().split('\t')
            a = int(row[0].strip()[4:])
            b = row[1].strip()
            ofh.writelines('\t'.join([str(a), b]) + '\n')


#~ ### FxnAss
filename = '/tmp/DrugCombination_2/FxnAss/STRINGTopo_EB_Allhuman.txt'
with open(filename) as ifh:
    with open(filename.replace('.txt', '_parsed.txt'), 'w') as ofh:
        for line in ifh:
            row = line.strip().split('\t')
            a = row[0].strip()
            b = int(row[1].strip().split('_')[0][4:])
            c = int(row[1].strip().split('_')[1][4:])
            d = row[2].strip()
            e = row[3].strip()
            f = row[4].strip()
            g = row[5].strip()
            if b <= c:
                ofh.writelines('\t'.join([a, str(b), str(c), d, e, f, g]) + '\n')
            else:
                ofh.writelines('\t'.join([a, str(c), str(b), d, e, f, g]) + '\n')


filename = '/tmp/DrugCombination_2/FxnAss/STRINGTopo_NSP_Allhuman.txt'
with open(filename) as ifh:
    with open(filename.replace('.txt', '_parsed.txt'), 'w') as ofh:
        for line in ifh:
            row = line.strip().split('\t')
            a = row[0].strip()
            b = int(row[1].strip().split('_')[0][4:])
            c = int(row[1].strip().split('_')[1][4:])
            d = row[2].strip()
            if b <= c:
                ofh.writelines('\t'.join([a, str(b), str(c), d]) + '\n')
            else:
                ofh.writelines('\t'.join([a, str(c), str(b), d]) + '\n')


filename = '/tmp/DrugCombination_2/FxnAss/STRINGTopo_Allhuman.txt'
with open(filename) as ifh:
    with open(filename.replace('.txt', '_parsed.txt'), 'w') as ofh:
        for line in ifh:
            row = line.strip().split('\t')
            a = row[0].strip()
            b = int(row[1].strip().split('_')[0][4:])
            c = int(row[1].strip().split('_')[1][4:])
            d = row[2].strip()
            e = row[3].strip()
            f = row[4].strip()
            g = row[5].strip()
            h = row[6].strip()
            i = row[7].strip()
            if b <= c:
                ofh.writelines('\t'.join([a, str(b), str(c), d, e, f, g, h, i]) + '\n')
            else:
                ofh.writelines('\t'.join([a, str(c), str(b), d, e, f, g, h, i]) + '\n')


#~ ### Drug
filename = '/tmp/DrugCombination_2/TargetAlias.txt'
with open(filename) as ifh:
        ifh.readline()
        with open(filename.replace('.txt', '_parsed.txt'), 'w') as ofh:
                for line in ifh:
                        row = line.strip().split('\t')
                        a = int(row[0].strip()[4:])
                        try:
                                b = int(row[1].strip()[4:]) if row[1].strip() != 'Not Available' else ''
                        except ValueError:
                                print line
                                print row
                                b = int(row[1].strip()[4:-1]) if row[1].strip() != 'Not Available' else ''
                        c = row[2].strip() if row[2].strip() != 'Not Available' else ''
                        d = int(row[3].strip()) if row[3].strip() != 'Not Available' else ''
                        e = row[4].strip() if row[4].strip() != 'Not Available' else ''
                        f = int(row[5].strip()) if row[5].strip() != 'Not Available' else ''
                        ofh.writelines('\t'.join([str(a), str(b), c, str(d), e, str(f)]) + '\n')


with open(filename) as ifh:
    ifh.readline()
    with open(filename.replace('.txt', '_parsed.txt'), 'w') as ofh:
        for line in ifh:
            row = line.strip().split('\t')
                try:
                    a = int(row[0].strip()[3:])
                    b = row[1].strip()
                    c = row[2].strip() if row[2].strip() != '-' else ''
                    d = row[3].strip()
                    e = row[4].strip()
                    f = row[5].strip()
                    g = row[6].strip()
                    h = row[7].strip()
                    i = row[8].strip() if row[8].strip() != '-' else ''
                    j = row[9].strip() if row[9].strip() != '-' else ''
                    ofh.writelines('\t'.join([str(a), b, c, d, e, f, g, h, i, j]) + '\n')
                except ValueError:
                    print line
                    print row


with open(filename) as ifh:
       ifh.readline()
       with open(filename.replace('.txt', '_parsed.txt'), 'w') as ofh:
               for line in ifh:
                       row = line.strip().split('\t')
                       a = row[0].strip()
                       b = int(row[1].strip().split('_')[0][4:])
                       c = int(row[1].strip().split('_')[1][4:])
                       d = int(row[2].strip().split('_')[0][3:])
