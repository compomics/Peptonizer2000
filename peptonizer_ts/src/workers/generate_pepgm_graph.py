import peptonizer

# The PSM input should be provided to the parser as a list of strings
psms = globals().get('input')

# print("Started parsing pout file from MS2Rescore...")
parsed_input = peptonizer.parse_ms2rescore_output(psms, 0.01)
# print(f"Input has been parsed successfully... --> size: {getsizeof(parsed_input)}")

# print("Started fetching Unipept taxon information...")
unipept_responses = peptonizer.fetch_unipept_taxon_information(
    parsed_input,
    "2",
    "species",
    "file_unipept_taxon_information_log"
)
# print("Successfully fetched Unipept taxon information...")

# print("Started weighing taxa...")
taxa_weights_df, _ = peptonizer.perform_taxa_weighing(
    unipept_responses,
    parsed_input,
    10,
    "species"
)
# print("Successfully weighed taxa...")

# print("Start creation of PepGM graph...")
pepgm_graph = peptonizer.generate_pepgm_graph(taxa_weights_df)
# print("Successfully created PepGM graph...")

pepgm_graph.to_graph_ml()
