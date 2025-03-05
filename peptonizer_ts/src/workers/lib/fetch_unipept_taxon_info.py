import peptonizer
import json

# The PSM input should be provided to the parser as a list of strings
pep_scores = globals().get('peptides_scores')
rank = globals().get('rank')
taxon_query = globals().get('taxon_query')

# First retrieve the taxonomic information from Unipept
unipept_responses = peptonizer.fetch_peptides_and_filter_taxa(
    list(pep_scores.keys()),
    taxon_query,
    rank,
    peptonizer.UnipeptCommunicator()
)

json.dumps(unipept_responses)
