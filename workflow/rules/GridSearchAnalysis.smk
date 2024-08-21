rule FindBestParameters:
     input: 
            ResultsDir + 'ClusterSortedTaxa.csv',
            expand(ResultsDir + 'Prior{Prior}/'+ 'PepGM_Results_a{alpha}_b{beta}_p{Prior}.png',alpha = AlphaRange, beta = BetaRange, Prior = prior)
     output: ResultsDir +'paramcheck.png', 
             ResultsDir + 'PeptonizerResults.csv'
             #ResultsDir +'PhyloTreeView.png'
     conda: 'envs/GridSearch.yml'
     log : ResultsDir +'paramcheck.log'
     params: Results = ResultsDir
     shell: 'python3 workflow/scripts/SearchAnalysisOutputFormatter.py --weights {input[0]} --resultsfolder {params.Results} --out {output[0]} &>> {log}'

