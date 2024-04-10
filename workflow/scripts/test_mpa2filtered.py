import requests
import json
import re
import os.path
from ete3 import NCBITaxa
import argparse
import logging
import time

import sys

ncbi = NCBITaxa()
url = 'http://api.unipept.ugent.be/mpa/pept2filtered'

def generatePostRequestChunks(peptides,TargetTaxa,chunksize=10,cutoff= 1000):
    '''
    Generates POST requests (json) queryong a chunk of peptides from petide list and target taxon
    :param peptides: list of peptides to query in Unipept
    :param TargetTaxa: list of one or more taxa to include in the Unipept query
    :param chunksize: number of peptides to be requested from Unipept
    :param cutoff: number of proteins a peptide is associated to above which said peptide will be removed from the query by Unipept. This enhances query speed.
    '''
    print('querying taxa ', TargetTaxa)
    AllTargetTaxa = []
    for Taxon in TargetTaxa:
        AllTargetTaxa.append(Taxon)
        AllTargetTaxa.extend(ncbi.get_descendant_taxa(Taxon, collapse_subspecies=False))
    
    
    Listofpeptides = [peptides[i:i + chunksize] for i in range(0, len(peptides), chunksize)]
    Listofrequests = [{"cutoff":cutoff, "peptides":chunk, "taxa":AllTargetTaxa} for chunk in Listofpeptides]

    return Listofrequests


requests_list = generatePostRequestChunks(['AEMNFNNPENGWFMDASVLFNNR','AQDLQLGFSYMF', 'HYQLSLSYSR',' YSTYSLLVPVNVGYK','METLYYK'],[1],chunksize=4,cutoff=1000)

for i in requests_list:
    request = requests.post(url,json.dumps(i),headers={'content-type':'application/json'}, timeout = 1800) 
    print(request.text)
