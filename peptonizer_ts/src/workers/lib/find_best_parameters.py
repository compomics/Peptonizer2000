import peptonizer
import networkx as nx
import pandas as pd



clustered_taxa_df = peptonizer.cluster_taxa_based_on_similarity(
    nx.read_graphml(args.full_graphml_path),
    pd.read_csv(args.taxa_weights_dataframe_file),
    args.similarity_threshold
)

