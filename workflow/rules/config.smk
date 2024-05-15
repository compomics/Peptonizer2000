configfile: 'config/config.yaml'


# global variables for each sample to import from tsv
#SAMPLES = pd.read_csv(config["samples"],index_col="sample",sep='\t')
#MGF_FILE = lambda wildcards: SAMPLES.at[wildcards.sample, 'mgf']
#Database = lambda wildcards: SAMPLES.at[wildcards.sample, 'Database']
#TaxaRank = lambda wildcards: SAMPLES.at[wildcards.sample, 'TaxaRank']

#Directories
DataDirectory = config['DataDir']
DatabaseDirectory = config ['DatabaseDir']
ExperimentDir = config['ResultsDir'] + config['ExperimentName'] 
ResultsDir = ExperimentDir+'/'+config['SampleName']+'/'
XTandemDir = ExperimentDir +'/'+config['SampleName']+'/XTandem/'
ExperimentDir = config['ResultsDir'] + config['ExperimentName']+'/' +config['SampleName']+'/'
MS2RescoreDir = config['ResultsDir'] + config['ExperimentName'] +'/'+config['SampleName']+'/MS2Rescore/'
ResourcesDir = config['ResourcesDir']
ResultsDirStrain = config['ResultsDir']

#Specific workflow settings
FilterSpectra = config['FilterSpectra']
AddCrapToDB = config['AddHostandCrapToDB']
ExperimentName = config['ExperimentName']
TaxaRank = config['TaxaRank']
InputSpectra = config['InputSpectra']
ReferenceDBName = config['ReferenceDBName']
SampleName = config['SampleName']
Pout = config['StartFromPout']
PoutFiles = config['PoutFiles']
StartFromUnipeptResponse = config['StartFromUnipeptResponse']
PreviousUnipeptQueryResults = config['PreviousUnipeptQueryResultsDir']

TaxaInPlot = config['TaxaInPlot']
#TaxaInProteinCount = config['TaxaInProteinCount']
#sourceDB = config['sourceDB']

AlphaRange = config['Alpha']
BetaRange = config['Beta']
prior = config['prior']
Regularize = config['Regularize']
print(Regularize, ' inside config smk')
# X!Tandem parameters
XTANDEM_DEFAULT = config['xtandem_default']
XTANDEM_FMME = config['xtandem_fmme']
XTANDEM_FMMEU = config['xtandem_fmmeu']
XTANDEM_PMMEP = config['xtandem_pmmep']
XTANDEM_PMMEM = config['xtandem_pmmem']
XTANDEM_PMMEU = config['xtandem_pmmeu']
# X!Tandem PTMs (comma separated if more than one)
XTANDEM_MODS_FIX = config['xtandem_mods_fixed'] 
XTANDEM_MODS_VAR = config['xtandem_mods_variable'] 
XTANDEM_MODS_VAR_NTERM = config['xtandem_mods_variable_nterm'] 
# additional parameters (dict: param name -> value)
XTANDEM_PARAS = config['xtandem_add_params']

#MS2rescore parameters
RescorePipeline = config['RescorePipeline']
RescoreFeatures = config['RescoreFeatures']
RunPercolator = config['RunPercolator']
FragModel = config['FragModel']
mods = config['Mods']
SpectrumTitlePattern = config['SpectrumTitlePattern']

#UnipeptQueryParameter
TaxaNumber = config['TaxaNumber']
targetTaxa = config['targetTaxa']
FDR = config['FDR']

#ParameterEvaluation Parameter
SimilarityThreshold = config['SimilarityThreshold']

#BiomassContributionPlotting
PosteriorThreshold = config['PosteriorThreshold']