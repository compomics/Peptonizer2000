import numpy as np
import pandas as pd
import os
import re 
import matplotlib
from matplotlib import pyplot as plt
from ete3 import NCBITaxa
from  scipy.stats import entropy
import rbo


#package settings
matplotlib.use('Agg')
ncbi = NCBITaxa()



def ComputeMetric(resultsfolder, output, weightsfile):

    """
    Compute a "goodness" metric for the parameters used in the grid search. 
    Generates a barplot of the metric and returns the best identified parameters.
    Return a list of the best parameter set in order [alpha,beta,gamma].

    :param resultsfolder: str, path to PepGM resultsfolder
    :param host: str, host to be excluded from the parameter checked taxa
    :param output: str, name of the output .png file 
    :param weightsfile: str, path to the .csv file with Unipept output
    
    """
    
    #predefine necessary lists
    Params = []
    Metrics = []
    Entropies = []
    #file with weights of taxids
    Weights = pd.read_csv(weightsfile,usecols=['HigherTaxa','scaled_weight'])
         
    for folders in os.listdir(resultsfolder):
        if os.path.isdir(resultsfolder+ '/' + folders):
            print(folders)
            for file in os.listdir(resultsfolder+ '/' + folders):
                print(file)
    
                if file.endswith('.csv'):
               
                    Results = pd.read_csv(resultsfolder+ '/' + folders +'/' + file,names = ['ID','score','type'])
                    TaxIDS = Results.loc[Results['type']=='taxon']
                    TaxIDS.loc[:,'score'] = pd.to_numeric(TaxIDS['score'],downcast = 'float')
                    TaxIDS.sort_values('score', ascending = False, inplace= True)
                    TaxIDS = TaxIDS.dropna()

    
                    #compute entropy of the posterior probability distribution
                    Entropy = entropy(TaxIDS['score'])
                    Entropies.append(Entropy)
                  
                    #compute the rank based similarity between the weight sorted taxa and the score sorted ID results
                    Matching = rbo.RankingSimilarity(Weights['HigherTaxa'].values,[int(i) for i in TaxIDS['ID'].values]).rbo()*(1/Entropy**2)
                    Metrics.append(1/Entropy**2)
    
                    #get the corresponding parameters
                    params = re.findall('\d+\.\d+',file)
                    Params.append([params[0],params[1],params[2]])


    zipLists = zip(Metrics,Params)
    SortedPairs = sorted(zipLists,reverse= True)
    
    tuples = zip(*SortedPairs)
    Metrics,Params = [list(tuple) for tuple in tuples]
    
    fig, ax = plt.subplots()
    fig.set_size_inches(30,15)
    
    ax.barh(list(range(len(Metrics[0:15]))),Metrics[0:15],color = 'mediumvioletred')
    plt.xticks(fontsize =35)
    plt.yticks(list(range(len(Metrics[0:15]))),Params[0:15],fontsize=25)

    plt.savefig(output)
    plt.close()

    return Params[0]




    




