import peptonizer
import json

# The PSM input should be provided to the parser as a list of strings
pep_scores = globals().get('peptides_scores')

# First retrieve the taxonomic information from Unipept
unipept_responses = peptonizer.fetch_unipept_taxon_information(
    list(pep_scores.keys()),
    "2",
    "species",
    "file_unipept_taxon_information_log"
)

json.dumps(unipept_responses)
