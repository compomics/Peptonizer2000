import argparse
import json
import numpy as np
from ete3 import NCBITaxa

import pandas as pd

ncbi = NCBITaxa()


def is_gzipped_file(filename):
    """ Check if input is gzipped."""
    with open(filename, 'rb') as f:
        head = f.read(2)
    return head == b'\x1f\x8b'


def init_argparser():
    """Init argument parser."""
    parser = argparse.ArgumentParser()

    parser.add_argument('--UnipeptResponseFile', type=str, required=True, help='path to Unipept response .json file')
    parser.add_argument('--NumberOfTaxa', type=int, required=True, help='number of taxa to include in the output')
    parser.add_argument('--out', type=str, required=True, help='path to csv out file')
    parser.add_argument('--TaxaWeightFile', type = str, required = False, help='path to taxaweight outputfile')
    parser.add_argument('--UnipeptPeptides', type=str, required=True, help='path to Unipept response .json file')
    parser.add_argument('--PeptidomeSize', type=str, required=False, help='path to proteome size per taxID file')
    parser.add_argument('--TaxaRank',type = str, required = False, default= 'species', help = 'taxonomic rank at which you want the peptonizer results resolved')

    args = parser.parse_args()

    return args


def GetPeptideCountPerTaxID(proteins_per_taxon):
    """
    Convert tab-separated taxon-protein counts to a dictionary.
    Parameters
    ----------
    proteins_per_taxon: str,
        Path to file that contains the counts.

    """
    protein_counts_per_taxid = {}

    with open(proteins_per_taxon, 'rb') as f:
        for line in f:
            # The file is encoded
            line_str = line.decode('utf-8').strip()
            # First column represents the TaxID, the second the count of peptides that are associated with that TaxID
            taxid, count = map(int, line_str.split('\t'))
            protein_counts_per_taxid[taxid] = int(count)

    return protein_counts_per_taxid

def GetLineageAtSpecifiedRank(taxid, TaxaRank):
    '''
    Returns the taxid of the specified rank in the lineage of any taxid given as argument
    -----
    taxid: int
         taxid to get the lineage of
    TaxaRank:
         rank you want to get
    '''
    taxids = set()
    for tax in taxid:
        RankTaxidDict = {}
        RankTaxidDict =  ncbi.get_rank(ncbi.get_lineage(tax))
        RankTaxidDict = {rank:taxiid for taxiid,rank in RankTaxidDict.items()}

        try:
            taxids.add(RankTaxidDict[str(TaxaRank)])
        except:
            taxids.add(tax)
        

    
    return list(taxids)


def WeightTaxa(UnipeptResponse, PeptScoreDict, MaxTax,*PeptidesPerTaxon, chunks=True, N=0, SelectRank = True, TaxaRank):
    """
    Weight inferred taxa based on their (1) degeneracy and (2) their proteome size.
    Parameters
    ----------
    UnipeptResponse: str
        Path to Unipept response json file
    PeptScoreDict: dict
        Dictionary that contains peptide to score & number of PSMs map
    MaxTax: int
        Maximum number of taxons to include in the graphical model
    PeptidesPerTaxon: str
        Path to the file that contains the size of the proteome per taxID (tab-separated)
    chunks: bool
        Allow memory-efficient reading of large json files
    N: int
        tbd

    Returns
    -------
    dataframe
        Top scoring taxa

    """
    with open(PeptScoreDict, 'r') as file:
        PeptScoreDictload = json.load(file)

    if chunks:
        with open(UnipeptResponse, 'r') as file:
            UnipeptDict = {"peptides": []}
            for line in file:
                try:
                    UnipeptDict["peptides"].extend(json.loads(line)["peptides"])
                except:
                    # TODO: Pieter fixes internal server error
                    # in the meantime, we work with the incomplete mapping
                    # UnipeptDict["peptides"] = [json.loads(line)["peptides"]]
                    continue

    else:
        with open(UnipeptResponse, 'r') as file:
            UnipeptDict = json.load(file)

    # Convert a JSON object into a Pandas DataFrame
    # record_path Parameter is used to specify the path to the nested list or dictionary that you want to normalize
    UnipeptFrame = pd.json_normalize(UnipeptDict, record_path=['peptides'])
    # Merge psm_score and number of psms
    UnipeptFrame = pd.concat([UnipeptFrame,
                              pd.json_normalize(UnipeptFrame['sequence'].map(PeptScoreDictload))], axis=1)
    # Score the degeneracy of a taxa, i.e.,
    # how conserved a peptide sequence is between taxa.
    #map all taxids in the list in the taxa column back to their taxid at species level
    UnipeptFrame['HigherTaxa'] = UnipeptFrame.apply(lambda row: GetLineageAtSpecifiedRank(row['taxa'],TaxaRank),axis = 1)

    # Divide the number of PSMs of a peptide by the number of taxa the peptide is associated with, exponentiated by 3
    UnipeptFrame['weight'] = UnipeptFrame['psms'].div([len(element)**3 for element in UnipeptFrame['HigherTaxa']])
    mask =[len(element)==1  for element in UnipeptFrame['HigherTaxa']]
    UniquePSMTaxa = set(i[0] for i in UnipeptFrame['HigherTaxa'][mask])
    UnipeptFrame = UnipeptFrame.explode('HigherTaxa', ignore_index=True)
    



    # Sum up the weights of a taxon and sort by weight
    UnipeptFrame['log_weight'] = np.log10(UnipeptFrame['weight']+1)
    TaxIDWeights = UnipeptFrame.groupby('HigherTaxa')['log_weight'].sum().reset_index()
    # Retrieve the proteome size per taxid as a dictionary
    # This file was previously prepared by filtering a generic accession 2 taxid mapping file
    # to swissprot (i.e., reviewed) proteins only
    

    #Peptidome size: optional to cinlude a weighting based on the size of the proteome, this didn't prove effective so disabled for now
    #PeptidomeSize = GetPeptideCountPerTaxID(PeptidesPerTaxon)
    # Map peptidome size and remove NAs
    #TaxIDWeights = TaxIDWeights[TaxIDWeights['taxa'].isin(PeptidomeSize.keys())].assign(
    #   proteome_size=lambda x: x['taxa'].map(PeptidomeSize))

    # Since large proteomes tend to have more detectable peptides,
    # we adjust the weight by dividing by the size of the proteome i.e.,
    # the number of proteins that are associated with a taxon
    TaxIDWeights["scaled_weight"] = TaxIDWeights["log_weight"] #/ (TaxIDWeights["proteome_size"]) ** N

   

    #retrieves the speciefied taxonomic rank taxid in the lineage of each of the species-level taxids returned by Unipept for both the UnipeptFrame 
    # and the TaxIdWeightFrame
    if SelectRank == True:
    #    TaxIDWeights['HigherTaxa'] = TaxIDWeights.apply(lambda row: GetLineageAtSpecifiedRank(row['taxa'],TaxaRank), axis = 1)
    #    UnipeptFrame['HigherTaxa'] = UnipeptFrame.apply(lambda row: GetLineageAtSpecifiedRank(row['taxa'],TaxaRank), axis = 1)
         HigherUniquePSMtaxids = UniquePSMTaxa#set([GetLineageAtSpecifiedRank(i,TaxaRank) for i in UniquePSMTaxa])
  
    
 
    
    #group the duplicate entries of higher up taxa and sum their weights    
    HigherTaxidWeights = TaxIDWeights.groupby('HigherTaxa')['scaled_weight'].sum().reset_index().sort_values(by=['scaled_weight'],
                                                                                          ascending=False)   
    #HigherTaxidWeights = TaxIDWeights
    HigherTaxidWeights['Unique'] = np.where(HigherTaxidWeights['HigherTaxa'].isin(HigherUniquePSMtaxids),True,False)

    try:
        HigherTaxidWeights = HigherTaxidWeights[HigherTaxidWeights.HigherTaxa!=1869227]
    except:
        pass

    if len(HigherTaxidWeights.HigherTaxa) < 50:
        return UnipeptFrame,HigherTaxidWeights
    else:
        TaxaToInclude = set(HigherTaxidWeights['HigherTaxa'][0:MaxTax])
        TaxaToInclude.update(HigherUniquePSMtaxids)
        UnipeptFrame[UnipeptFrame['HigherTaxa'].isin(TaxaToInclude)]
        return UnipeptFrame[UnipeptFrame['HigherTaxa'].isin(TaxaToInclude)],HigherTaxidWeights

args = init_argparser()
DF,Weights = WeightTaxa(args.UnipeptResponseFile, args.UnipeptPeptides, args.NumberOfTaxa, TaxaRank = args.TaxaRank)
DF.to_csv(args.out)
Weights.to_csv(args.TaxaWeightFile)

#if __name__ == '__main__':
    #args = init_argparser()
    #DF = WeightTaxa(args.UnipeptResponseFile, args.UnipeptPeptides, args.NumberOfTaxa, args.PeptidomeSize, TaxaRank = args.TaxaRank)
    #DF.to_csv(args.out)

   
    #DF = WeightTaxa('/home/tanja/Peptonizer2000/Peptonizer2000/results/CAMPI1_SIHUMIx_allbacteria/CAMPI_SIHUMIx/UnipeptResponse.json','/home/tanja/Peptonizer2000/Peptonizer2000/results/CAMPI1_SIHUMIx_allbacteria/CAMPI_SIHUMIx/UnipeptPeptides.json',
                    #50,'/home/tanja/Peptonizer2000/Peptonizer2000/resources/taxa_peptidome_size.tsv')
    #DF.to_csv('/home/tanja/Peptonizer2000/Peptonizer2000/results/CAMPI1_SIHUMIx_allbacteria/CAMPI_SIHUMIx/GraphDataframe_test.csv')
    