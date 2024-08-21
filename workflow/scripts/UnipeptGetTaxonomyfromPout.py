import requests
import json
import re
import os.path
from ete3 import NCBITaxa
import argparse
import logging
import time
import asyncio
import aiohttp
import aiofiles
import sys

ncbi = NCBITaxa()

parser = argparse.ArgumentParser()

parser.add_argument('--UnipeptResponseFile', type = str, required = True, help = 'path to Unipept response .json file')
parser.add_argument('--TaxonomyQuery', required = True, help = 'taxa to query in Unipept. If querying all taxa, put [1]') 
parser.add_argument('--PeptidesAndScores', type = str, required = True, help = 'path to tsv file with peptides and scores')
parser.add_argument('--pep_out', type = str, required = True, help = 'path to csv out file')
parser.add_argument('--logfile', type = str, required = True, help ='path to logfile of failed unipeptquery attempts')
#parser.add_argument('--UnipeptRequestLim', type = int, required = False, default = 3 , help = ' number of requests allowed to be executed in parallel on the UNipept server')


args = parser.parse_args()


def parse_tsv_to_peptide_dict(tsv_file):
    '''
    Parses a .tsv file for peptides and scores.
    
    :param tsv_file: str, path to the tsv file
    :return: dict, peptides: [highest_score, #occurrences]
    '''
    
    peptide_dict = {}
    
    with open(tsv_file, "r") as f:
        next(f)  # Skip header
        for line in f:
            if line.rstrip() == "":
                continue
            
            splitted_line = line.rstrip().split("\t")
            
            if len(splitted_line) < 2:
                continue  # Skip lines that don't have at least two columns
            
            peptide, score = splitted_line[0], float(splitted_line[1])
            
            # Replace score of 1 with 0.99
            if score == 1:
                score = 0.99
            
            if peptide not in peptide_dict:
                peptide_dict[peptide] = [score, 1]
            else:
                # Update the highest score
                peptide_dict[peptide][0] = max(peptide_dict[peptide][0], score)
                # Increment occurrence count
                peptide_dict[peptide][1] += 1
    
    return peptide_dict


def PepListNoMissedCleavages(peptide):
    """
    Takes a peptide and cleaves it into Unipept format (0 missed cleavages,
    cleaves after K or R except followed by P)
    :param peptides_in: list of peptides
    :return: list of peptides in Unipept format
    """
    
    peptides = list()
    trypsin = lambda peptide: re.sub(r'(?<=[RK])(?=[^P])', '\n', peptide, re.DOTALL).split()
    peptides += trypsin(peptide)

    return peptides


def generatePostRequestChunks(peptides,TargetTaxa,chunksize=10,cutoff= 1000):
    '''
    Generates POST requests (json) queryong a chunk of peptides from petide list and target taxon
    :param peptides: list of peptides to query in Unipept
    :param TargetTaxa: list of one or more taxa to include in the Unipept query
    :param chunksize: number of peptides to be requested from Unipept
    :param cutoff: number of proteins a peptide is associated to above which said peptide will be removed from the query by Unipept. This enhances query speed.
    '''
    print('querying taxa ', TargetTaxa)
    AllTargetTaxa = set()
    for Taxon in TargetTaxa:
        AllTargetTaxa.add(Taxon)
        AllTargetTaxa.update(ncbi.get_descendant_taxa(Taxon, collapse_subspecies=False))
        AllTargetTaxa.update(ncbi.get_descendant_taxa(Taxon, collapse_subspecies=True))

    
    
    Listofpeptides = [peptides[i:i + chunksize] for i in range(0, len(peptides), chunksize)]
    Listofrequests = [{"cutoff":cutoff, "peptides":chunk, "taxa":list(AllTargetTaxa)} for chunk in Listofpeptides]

    return Listofrequests



def PostInfoFromUnipeptChunks(request_json, out_file, failed_requests_file):
    """
    Send all requests, get for each peptide the phylum, family, genus and collection of EC-numbers
    :param request_list: list of Get Requests
    :param result_file: csv file with Unipept info (phylum, family, genus and collection of EC-numbers)
    :return: None
    """

    logging.basicConfig(filename= failed_requests_file, level=logging.INFO)
    
    url = 'http://api.unipept.ugent.be/mpa/pept2filtered'
    print('now querying Unipept in '+str(len(request_json))+' chunks')
    
    
    #try Unipept query in chunks
    failed_requests = {}
    count = 0
    for chunk in request_json:
        start = time.time()
        try:
            
            request = requests.post(url,json.dumps(chunk),headers={'content-type':'application/json'}, timeout = 1800) 
            request.raise_for_status()

            with open(out_file, 'a') as f_out:
                print(request.text,file=f_out)
            end = time.time()
            print('sucessfully queried a chunk, number '+str(count) +' out of '+str(len(request_json))+'\n time taken' + str(end-start))
            count += 1

        
        except requests.exceptions.RequestException as e:

            logging.error(f'Request {url} failed with error {e}')
            failed_requests[json.dumps(chunk)] = e

    #retry failed requests        
    for chunk,error in failed_requests.items():
        try: 
            request = requests.post(url,chunk,headers={'content-type':'application/json'}, timeout = 3600) 
            request.raise_for_status()
            with open(out_file, 'a') as f_out:
                print(request.text,file=f_out)
        
        except requests.exceptions.RequestException as e:
            logging.error(f'Retry request to {url} failed with error: {e}')


async def fetch_data(session, url, json_input, out_file, failed_requests, i, total_chunks):
    try:
        async with session.post(url, json=json_input, timeout=8600, headers={'content-type': 'application/json'}) as response:
            response.raise_for_status()
            result = await response.text()
            async with aiofiles.open(out_file, 'a') as f_out:
                await f_out.write(result + '\n')
            print(f'Successfully queried chunk {i} out of {total_chunks}')
            return True
    except (aiohttp.ClientError,aiohttp.ClientConnectionError, asyncio.TimeoutError) as e:
        print(f'Failed to query chunk {i}')
        logging.error(f'Request to {url} failed with error: {e}')
        failed_requests.append((json_input, e))
        return False

async def limited_gather(semaphore, session, url, chunks, out_file, failed_requests, i, total_chunks):
    async with semaphore:
        return await fetch_data(session, url, chunks, out_file, failed_requests, i, total_chunks)

async def main(request, out_file, failed_requests_file):
    logging.basicConfig(filename=failed_requests_file, level=logging.INFO)
    url = 'http://api.unipept.ugent.be/mpa/pept2filtered'
    print(f'Now querying Unipept in {len(request)} chunks')

    semaphore = asyncio.Semaphore(3)  # Limit the number of concurrent requests to three

    total_chunks = len(request)
    failed_requests = []

    async with aiohttp.ClientSession() as session:
        limited_gather_tasks = []
        for i, chunk in enumerate(request):
            task = limited_gather(semaphore, session, url, chunk, out_file, failed_requests, i, total_chunks)
            limited_gather_tasks.append(task)

        await asyncio.gather(*limited_gather_tasks)

        # Retry failed requests with exponential backoff
        for retry in range(3):  # Retry a maximum of 3 times
            retry_tasks = []
            for chunk, e in failed_requests:
                logging.error(f'Request {chunk} try 2 to {url} failed with error: {e}')
                task = limited_gather(semaphore, session, url, chunk, out_file, failed_requests, i, total_chunks)
                retry_tasks.append(task)

            await asyncio.gather(*retry_tasks)
            await asyncio.sleep(2 ** retry)  # Exponential backoff


pep_score_psm = parse_tsv_to_peptide_dict(args.PeptidesAndScores)
UnipeptPeptides = dict()
for peptide in pep_score_psm.keys():
    FullyTrypticPeptides = PepListNoMissedCleavages(peptide)
    for pep in FullyTrypticPeptides:
        UnipeptPeptides[pep] ={'score':pep_score_psm[peptide][0], 'psms':pep_score_psm[peptide][1]} 
with open(args.pep_out, 'a') as f_out:
    f_out.write( json.dumps(UnipeptPeptides))

#get and save Info from Unipept if the response file doesn't exist yet
request = generatePostRequestChunks(list(UnipeptPeptides.keys()),[int(item) for item in args.TaxonomyQuery.split(',')])    
#save = PostInfoFromUnipeptChunks(request,args.UnipeptResponseFile,args.logfile)
save = asyncio.run(main(request, args.UnipeptResponseFile, args.logfile))


#if __name__=='__main__':
#    pout_file = '/home/tholstei/repos/PepGM_all/PepGM/results/Xtandem_rescore_test/PXD018594_Sars_CoV_2/MS2Rescore/rescored_searchengine_ms2pip_rt_features.pout'
#    pep_score_psm = Poutparser(pout_file,0.05,'')


 #   UnipeptPeptides = dict()
 #   for peptide in pep_score_psm.keys():
 #       FullyTrypticPeptides = PepListNoMissedCleavages(peptide)
 #       for pep in FullyTrypticPeptides:
 #           UnipeptPeptides[pep] = pep_score_psm[peptide]

  #  save = asyncio.run(main(request, args.UnipeptResponse, args.logfile))


#    out = PostInfoFromUnipept(request,'test_unipept.json')
