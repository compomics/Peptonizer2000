from GridSearchAnalysis import *
#from PhyloTreeView import *
import pandas as pd
import os
import argparse 
from ete3 import NCBITaxa

parser = argparse.ArgumentParser(description = 'Downstream analysis of the grid search for the PepGM algorithm')
parser.add_argument('--resultsfolder', required = True, help = 'folder with the results from the beliefpropagation in csv format')
parser.add_argument('--out', required=True, help ='output png file')
parser.add_argument('--reference_db',required = True, help = 'name of the reference datasbe used to include in the filepath string')
parser.add_argument('--weights', required = True, help = 'path to the file of weighted taxids')

args = parser.parse_args()

ncbi = NCBITaxa()

matplotlib.use('Agg')


def SaveReducedCSV(filepath,output):
    '''
    Saves CSV of the 20 best results to a file in the main results folder
    :param filepath: str, path of cvs to be reduced
    :param output: str,output path of csv
    '''
                   
    Results = pd.read_csv(filepath,names = ['ID','score','type'])
    TaxIDS = Results.loc[Results['type']=='taxon']
    TaxIDS.loc[:,'score'] = pd.to_numeric(TaxIDS['score'],downcast = 'float')
    TaxIDS = TaxIDS.sort_values('score', ascending = False)
   
    TaxIDS.to_csv(output)
           
#SaveReducedCSV('/home/tanja/Peptonizer2000/Peptonizer2000/results/S07_IGC_reg/CAMPI_SIHUMIx/Prior0.5/GUT_DB1IGC_PepGM_Results_a0.99_b0.8_p0.5.csv','/home/tanja/Peptonizer2000/Peptonizer2000/results/S07_IGC_reg/CAMPI_SIHUMIx/PeptonizerResults.csv')

def MoveBestResultsPlot(filepath,out):
    '''
    Copies barplot of best identifiedparameter set to main resultsfolder
    :param filepath: list, the three grid search parameters identified as best for the sample at hand
    :param out: output path where the file will be copied

    '''
    os.system('cp '+filepath +' '+out)


#analyse the PepGM grid search with an empirical metric, retrieve best parameters
Parameters = ComputeMetric(args.resultsfolder, args.out,args.weights)
Resultsfile = args.resultsfolder+'/Prior'+str(Parameters[2]) +'/'+ args.reference_db + '_PepGM_Results_a'+str(Parameters[0])+'_b'+str(Parameters[1])+'_p'+str(Parameters[2])

#save the reduced results csv
SaveReducedCSV(Resultsfile+'.csv', args.resultsfolder +'/PeptonizerResults.csv')
MoveBestResultsPlot(Resultsfile+'.png',args.resultsfolder+'/PeptonizerResults.png')

#save a phylogenetic tree view of the PepGm results
#CreatePhyloTreeView(Resultsfile+'.csv',args.resultsfolder + '/PhyloTreeView.png')