  rule UnipeptQuery:
    input: 
          PeptidesAndScores
    params:
          targetTaxa = targetTaxa,
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