import networkx as nx
import numpy as np
import pandas as pd

CSVPath = '/home/tanja/Peptonizer2000/Peptonizer2000/results/CAMPI1_SIHUMIx_allbacteria_s03_new_rescore/CAMPI_SIHUMIx/GraphDataframe.csv'
peptonizer_results_file = '/home/tanja/Peptonizer2000/Peptonizer2000/results/CAMPI1_SIHUMIx_allbacteria_s03_new_rescore/CAMPI_SIHUMIx/PeptonizerResults.csv'
out = '/home/tanja/Peptonizer2000/Peptonizer2000/results/CAMPI1_SIHUMIx_allbacteria_s03_new_rescore/CAMPI_SIHUMIx/PeptonizerResults_greedy.csv'

def greedy_peptide_aggregation(CSVPath, peptonizer_results_file,out):
    UnipeptResponse = pd.read_csv(CSVPath)
    newGraph = nx.from_pandas_edgelist(UnipeptResponse, 'sequence', 'HigherTaxa')
    PeptideAttributes = UnipeptResponse.apply(lambda row: (row["sequence"], {'InitialBelief_0': row["score"],'InitialBelief_1': 1-row['score'],'category':'peptide'}) ,axis=1)
    TaxaAttributes = UnipeptResponse.apply(lambda row: (row["HigherTaxa"], {'category':'taxon'}) ,axis=1)
    G = nx.Graph()
    G.add_edges_from(newGraph.edges)
    G.add_nodes_from(PeptideAttributes)
    G.add_nodes_from(TaxaAttributes)    
    #load peptonizer results
    peptonizer_results = pd.read_csv(peptonizer_results_file)
    #get list of 'ID' column entries in peptonizer results
    taxon_ids = peptonizer_results['ID'].tolist()
    #for each node in the list of taxa in peptonizer results, get peptide nodes connected to it, count them, and remove them from graph
    taxon_peptide_count = []
    for taxon_id in taxon_ids:
        print(taxon_id)
        taxon_peptide_nodes = [node for node in G.neighbors(taxon_id) if G.nodes[node]['category'] == 'peptide']
        taxon_peptide_count.append(len(taxon_peptide_nodes))
        G.remove_nodes_from(taxon_peptide_nodes)

    peptonizer_results['peptide_count'] = taxon_peptide_count
    peptonizer_results.to_csv(out, index=False)

greedy_peptide_aggregation(CSVPath, peptonizer_results_file,out)

