###rule for parsing percolator file

# check if spectrum should be filtered or not
def PoutToUse(condition1,condition2):
    if condition1 and condition2:
        return expand('{poutfiles}.filtered', poutfiles = PoutFiles) 
    if condition1 and not condition2:
        return expand('{poutfiles}', poutfiles = PoutFiles) 
    if not condition1 and condition2:
        return expand(ExperimentDir+'{spectrum_name}/ms2rescore/rescored/rescored.filtered.psms.tsv',spectrum_name = SpectraNames)
    if not condition1 and not condition2:
        return expand(ExperimentDir+'{spectrum_name}/ms2rescore/rescored/rescored.psms.tsv',spectrum_name = SpectraNames)


def FilteredOutputPout(condition1):
    if condition1:
        return '{poutfiles}.filtered'
    else:
        return ExperimentDir+'{spectrum_name}/ms2rescore/rescored/rescored.filtered.psms.tsv'

def InputPout(condition1):
    if condition1:
        return '{poutfiles}'
    else:
        return ExperimentDir+'{spectrum_name}/ms2rescore/rescored/rescored.psms.tsv'

InputPoutFile = PoutToUse(Pout,FilterSpectra)
FilteredInputPout = InputPout(Pout)
OutputPoutFile = FilteredOutputPout(Pout)

rule FilterPout:
    input: 
          FilteredInputPout
    output:
          OutputPoutFile
    shell: "python3 workflow/scripts/FilterPout.py --Pout {input} --out {output}"
    
rule UnipeptQuery:
    input: 
          InputPoutFile
    params:
          targetTaxa = targetTaxa,
          FDR = FDR
    log: ResultsDir + 'UnipeptResponse.log'
    output: 
          ResultsDir + 'UnipeptResponse.json',
          ResultsDir + 'UnipeptPeptides.json'
    conda: 'envs/Unipeptquery.yml'   
    shell: "python3 workflow/scripts/UnipeptGetTaxonomyfromPout.py --UnipeptResponseFile {output[0]} --pep_out {output[1]} --TaxonomyQuery {params.targetTaxa} --FDR {params.FDR} --PoutFile {input} --logfile {log}" 


def StartFromUnipept(condition):
    if condition:
        return [PreviousUnipeptQueryResults + 'UnipeptResponse.json', PreviousUnipeptQueryResults + 'UnipeptPeptides.json']
    else:
        return [ResultsDir + 'UnipeptResponse.json', ResultsDir + 'UnipeptPeptides.json']



rule ParseToUnipeptCSV:
    input: 
          StartFromUnipept(StartFromUnipeptResponse)
          
          
    params: 
      NumberofTaxa = TaxaNumber,
      TaxaRank = TaxaRank
           
    log: ResultsDir + 'ParsetoCSV.log'
    output: 
            ResultsDir + 'GraphDataframe.csv',
            ResultsDir +'TaxaWeights.csv'
    conda: 'envs/graphenv.yml' 
    shell: "python3 workflow/scripts/WeightTaxa.py --UnipeptResponseFile {input[0]} --UnipeptPeptides {input[1]} --out {output[0]} --TaxaWeightFile {output[1]} --NumberOfTaxa {params.NumberofTaxa} --TaxaRank {params.TaxaRank}" 