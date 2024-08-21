configfile: 'config/config.yaml'

#Directories
DataDirectory = config['DataDir']
ExperimentDir = config['ResultsDir'] + config['ExperimentName'] 
ResultsDir = ExperimentDir+'/'+config['SampleName']+'/'
ExperimentDir = config['ResultsDir'] + config['ExperimentName']+'/' +config['SampleName']+'/'
ResourcesDir = config['ResourcesDir']
ResultsDirStrain = config['ResultsDir']

#Specific workflow settings
ExperimentName = config['ExperimentName']
TaxaRank = config['TaxaRank']
PeptidesAndScores = config['PeptidesAndScores']
SampleName = config['SampleName']
StartFromUnipeptResponse = config['StartFromUnipeptResponse']
PreviousUnipeptQueryResults = config['PreviousUnipeptQueryResultsDir']

TaxaInPlot = config['TaxaInPlot']

AlphaRange = config['Alpha']
BetaRange = config['Beta']
prior = config['prior']
Regularize = config['Regularize']
print(Regularize, ' inside config smk')


#UnipeptQueryParameter
TaxaNumber = config['TaxaNumber']
targetTaxa = config['targetTaxa']
FDR = config['FDR']

#ParameterEvaluation Parameter
SimilarityThreshold = config['SimilarityThreshold']

#BiomassContributionPlotting
PosteriorThreshold = config['PosteriorThreshold']