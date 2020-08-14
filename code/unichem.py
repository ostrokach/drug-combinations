# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 16:24:38 2013

@author: alexey


Query `unichem` to convert different chemical ids.
"""
import json
import urllib2
import numpy as np


who_drug_name = 'Metyrapone'
who_drug_name = who_drug_name.lower()
query_tuple = (who_drug_name, '12', '')

pubchem_cid = '152951'
query_tuple = (pubchem_cid, '22', '1')

try:
    request_string = 'https://www.ebi.ac.uk/unichem/rest/mapping/22/1'
    page = urllib2.urlopen(request_string).read()
except urllib2.HTTPError as err:
    print err
    print pubchem_cid, 'not found!'
    page =  "[{\"assignment\":\"0\",\"src_compound_id\":\"\"}]"

page2 = json.loads(page)
print page
print page2
if np.size(page2) > 0:
    for p in page2:
        print str(p.values())

        page2[0]['src_compound_id']

