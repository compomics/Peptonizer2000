from typing import List, Dict

import numpy as np
import pandas as pd

from .taxon_manager import TaxonManager


def get_lineage_at_specified_rank(tax_ids: List[int], taxa_rank: str) -> List[int]:
    """
    Returns the taxon ID of the specified rank in the lineage for all taxa ID given as argument.

    For example, given a taxon ID at strain level and "species" as value for the taxa_rank argument, this function
    will return the taxon ID at species level for the input taxon ID.

    Parameters
    -----
    tax_ids: [int]
         List of taxon_ids to get the lineage of
    taxa_rank:
         Rank at which you want to pin the taxa
    """

    # Get the full lineage for all given taxon IDs from Unipept
    lineages = TaxonManager.get_lineages_for_taxa(tax_ids)

    # Get the index of the NCBI rank that we're interested in. This index is required to extract the taxon IDs from the
    # correct place in the lineage.
    rank_idx = TaxonManager.NCBI_RANKS.index(taxa_rank)

    return [lineages[tax][rank_idx] for tax in tax_ids]


def weighted_random_sample(peptide_taxa: Dict[str, List[int]], n: int) -> Dict[str, List[int]]:
    """
    Randomly select n pairs from the provided peptide_taxa dictionary. The chance for each peptide to be selected,
    depends on the amount of taxa that are associated to this peptide. The more taxa, the lower the information
    content of this peptide, and the lower the chance of being selected.

    :param peptide_taxa: A dictionary mapping peptide sequences onto the taxon IDs that they're associated with.
    :param n: The amount of peptide-taxa pairs that should be selected from the provided peptide_taxa dictionary.
    :return: A new peptide_taxa dictionary object containing only the n selected pairs.
    """

    # Calculate weights based on the length of the taxa array for each peptide
    weights = np.array([1 / len(taxa) if taxa else 0 for taxa in peptide_taxa.values()])

    # Normalize weights
    total_weight = np.sum(weights)
    if total_weight == 0:
        raise ValueError("All objects have zero weight, cannot sample.")

    normalized_weights = weights / total_weight

    # Perform weighted sampling without replacement (so no items are selected more than once)
    peptides = list(peptide_taxa.keys())
    sampled_indices = np.random.choice(len(peptides), size=min(np.count_nonzero(normalized_weights), n), replace=False, p=normalized_weights)

    output = dict()

    # Retrieve sampled objects
    for sampled_index in sampled_indices:
        sampled_peptide = peptides[sampled_index]
        output[sampled_peptide] = peptide_taxa[sampled_peptide]

    return output


def normalize_taxa(peptide_taxa: Dict[str, List[int]], taxa_rank: str):
    """
    Map all taxon IDs that are found (as keys) in the given peptide_taxa dictionary to the parent or child at taxa_rank.

    :param peptide_taxa: A dictionary that maps peptide sequences onto a list of taxon IDs that are associated to this
        peptide.
    :param taxa_rank: The NCBI taxon rank at which the taxon IDs should be mapped. Must be a valid NCBI rank that's
        supported by Unipept.
    :return: The same dictionary that was provided as input, but then with the modified taxa. The returned object is the
        same as the input object.
    """
    # Map all taxa onto the rank specified by the user
    for (peptide, taxa) in peptide_taxa.items():
        peptide_taxa[peptide] = list(set(get_lineage_at_specified_rank(taxa, taxa_rank)))
    return peptide_taxa


def perform_taxa_weighing(
    peptide_taxa: Dict[str, List[int]],
    pep_scores: Dict[str, float],
    pep_psm_counts: Dict[str, int],
    max_taxa,
    taxa_rank="species"
):
    """
    Weight inferred taxa based on their (1) degeneracy and (2) their proteome size.

    Parameters
    ----------
    peptide_taxa: List[any]
        Peptide counts that have already been processed by Unipept before.
    pep_scores: Dict[str, Dict[str, float | int]]
        Dictionary that maps each peptide string onto an object containing the score associated to this peptide and the
        psm count.
    max_taxa: int
        Maximum number of taxa to include in the final graphical model.
    taxa_rank: str
        NCBI rank at which the Peptonizer analysis should be performed.

    Returns
    -------
    dataframe
        Top scoring taxa

    """
    print("Started mapping all taxon ids to the specified rank...")
    peptide_taxa = normalize_taxa(peptide_taxa, taxa_rank)
    peptide_taxa = weighted_random_sample(peptide_taxa, 10000)

    print(f"Using {len(peptide_taxa)} sequences as input...")

    print("Unipept responses:")
    print(peptide_taxa)

    # Convert a JSON object into a Pandas DataFrame
    # record_path Parameter is used to specify the path to the nested list or dictionary that you want to normalize
    print("Normalizing peptides and converting to dataframe...")
    unipept_frame = pd.DataFrame(list(peptide_taxa.items()), columns=['sequence', 'taxa'])

    scores = unipept_frame["sequence"].map(pep_scores)
    scores.name = "score"

    psms = unipept_frame["sequence"].map(pep_psm_counts)
    psms.name = "psms"

    # Merge psm_score and number of psms
    unipept_frame = pd.concat(
        [
            unipept_frame,
            scores,
            psms
        ],
        axis=1,
    )

    # Score the degeneracy of a taxa, i.e.,
    # how conserved a peptide sequence is between taxa.
    # map all taxids in the list in the taxa column back to their taxid at species level (or the rank specified by the user)
    # Right now, HigherTaxa is simply a copy of taxa. This step still needs to be optimized.
    unipept_frame["HigherTaxa"] = unipept_frame.apply(
        lambda row: row["taxa"], axis=1
    )

    # unipept_frame.to_csv("higher_taxa_step.csv")

    # Divide the number of PSMs of a peptide by the number of taxa the peptide is associated with, exponentiated by 3
    print("Started dividing the number of PSMS of a peptide by the number the peptide is associated with...")
    unipept_frame["weight"] = unipept_frame["psms"].div(
        [len(element) ** 3 for element in unipept_frame["HigherTaxa"]]
    )

    # unipept_frame.to_csv("weight_step.csv")

    mask = [len(element) == 1 for element in unipept_frame["HigherTaxa"]]
    unique_psm_taxa = set(i[0] for i in unipept_frame["HigherTaxa"][mask])
    unipept_frame = unipept_frame.explode("HigherTaxa", ignore_index=True)

    # Sum up the weights of a taxon and sort by weight
    print("Started summing the weights of a taxon and sorting them by weight...")
    unipept_frame["log_weight"] = np.log10(unipept_frame["weight"] + 1)
    tax_id_weights = unipept_frame.groupby("HigherTaxa")["log_weight"].sum().reset_index()

    # Since large proteomes tend to have more detectable peptides,
    # we adjust the weight by dividing by the size of the proteome i.e.,
    # the number of proteins that are associated with a taxon
    tax_id_weights["scaled_weight"] = tax_id_weights[
        "log_weight"
    ]  # / (TaxIDWeights["proteome_size"]) ** N

    # Retrieves the specified taxonomic rank taxid in the lineage of each of the species-level taxids returned by
    # Unipept for both the UnipeptFrame and the TaxIdWeightFrame
    higher_unique_psm_taxids = unique_psm_taxa  # set([GetLineageAtSpecifiedRank(i,TaxaRank) for i in UniquePSMTaxa])

    # group the duplicate entries of higher up taxa and sum their weights
    print("Started grouping duplicate entries of taxa situated higher up and sum their weights...")
    higher_taxid_weights = (
        tax_id_weights.groupby("HigherTaxa")["scaled_weight"]
        .sum()
        .reset_index()
        .sort_values(by=["scaled_weight"], ascending=False)
    )
    # HigherTaxidWeights = TaxIDWeights
    higher_taxid_weights["Unique"] = np.where(
        higher_taxid_weights["HigherTaxa"].isin(higher_unique_psm_taxids), True, False
    )

    try:
        higher_taxid_weights = higher_taxid_weights[
            higher_taxid_weights.HigherTaxa != 1869227
        ]
    except:
        pass

    unipept_frame = unipept_frame.drop('taxa', axis=1)

    if len(higher_taxid_weights.HigherTaxa) < 50:
        return unipept_frame, higher_taxid_weights
    else:
        taxa_to_include = set(higher_taxid_weights["HigherTaxa"][0:max_taxa])
        taxa_to_include.update(higher_unique_psm_taxids)
        return (
            unipept_frame[unipept_frame["HigherTaxa"].isin(taxa_to_include)],
            higher_taxid_weights,
        )
