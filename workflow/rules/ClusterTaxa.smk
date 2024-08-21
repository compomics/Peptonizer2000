rule ClusterTaxa:
    input: 
        ResultsDir + 'PepGM_graph.graphml',
        ResultsDir +'TaxaWeights.csv'
    params:
        SimThreshold = 0.9
    log: ResultsDir + 'ClusterTaxa.log'
    output: 
        ResultsDir + 'TaxaSimilarityMatrix.csv',
        ResultsDir + 'ClusterSortedTaxa.csv'
    conda: 'envs/graphenv.yml'   
    shell: "python3 workflow/scripts/taxa_clustering.py --outSimilarities {output[0]} --out {output[1]} --SimilarityThreshold {params.SimThreshold}  --TaxaWeightFile {input[1]} --GraphIN {input[0]}" 

