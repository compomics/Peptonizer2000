from typing import Dict, List, Any

from .unipept_communicator import UnipeptCommunicator

def fetch_peptides_and_filter_taxa(
    peptides: List[str],
    taxonomy_query: str,
    rank: str
) -> Dict[str, List[int]]:
    # First we retrieve all taxa associated with the given taxa
    unipept_communicator = UnipeptCommunicator()
    peptides_taxa = unipept_communicator.get_taxa_for_peptides(peptides)

    # Then, we make sure to filter the taxa and only keep those that are associated to the taxa of interest indicated by
    # the user.

    # Retrieve all (in)direct children of the filter taxa provided by the user
    taxa_filter = set(unipept_communicator.get_descendants_for_taxa([int(item) for item in taxonomy_query.split(",")], rank))

    # Compute the intersection of the taxa that should be retained and the original list of taxa
    for (peptide, taxa) in peptides_taxa.items():
        peptides_taxa[peptide] = list(set(taxa) & taxa_filter)

    return peptides_taxa


