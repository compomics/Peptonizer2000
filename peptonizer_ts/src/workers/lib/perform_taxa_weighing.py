import peptonizer

# The PSM input should be provided to the parser as a list of strings
pep_scores = globals().get('peptides_scores')
pep_counts = globals().get('peptides_counts')

# First retrieve the taxonomic information from Unipept
unipept_responses = peptonizer.fetch_unipept_taxon_information(
    list(pep_scores.keys()),
    "2",
    "species",
    "file_unipept_taxon_information_log"
)

# Infer the taxa weights for these peptide sequences
sequence_scores_df, taxa_weights_df = peptonizer.perform_taxa_weighing(
    unipept_responses,
    pep_scores,
    pep_counts,
    10,
    "species"
)

# Return a CSV-representation of the taxa_weights dataframe
output = [sequence_scores_df.to_csv(), taxa_weights_df.to_csv()]
output