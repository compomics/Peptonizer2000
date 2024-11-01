import logging
import requests

from typing import List, Dict
from sys import getsizeof

from .request_manager import RequestManager

class TaxonManager:
    UNIPEPT_URL = "http://api.unipept.ugent.be"
    UNIPEPT_TAXONOMY_ENDPOINT = "/api/v2/taxonomy.json"

    TAXONOMY_ENDPOINT_BATCH_SIZE = 100

    NCBI_RANKS = [
        "superkingdom",
        "kingdom",
        "subkingdom",
        "superphylum",
        "phylum",
        "subphylum",
        "superclass",
        "class",
        "subclass",
        "superorder",
        "order",
        "suborder",
        "infraorder",
        "superfamily",
        "family",
        "subfamily",
        "tribe",
        "subtribe",
        "genus",
        "subgenus",
        "species_group",
        "species_subgroup",
        "species",
        "subspecies",
        "strain",
        "varietas",
        "forma"
    ]

    # Static cache to store retrieved lineages
    lineage_cache = {}

    @staticmethod
    def get_descendants_for_taxa(target_taxa: List[int], descendants_rank: str) -> List[int]:
        """
        Returns a list of all taxon IDs that are a descendant of one of the given taxa (in the target_taxa argument).

        Parameters
        ----------
        target_taxa: List[int],
            A list of taxon IDs for which all descendants at a specific NCBI rank (and lower) should be retrieved
        descendants_rank: str,
            The maximum rank that each of the descendants should have in the NCBI taxonomy. Thus, all descendants that
            are defined at this rank, or deeper, are reported.
        """
        url = TaxonManager.UNIPEPT_URL + TaxonManager.UNIPEPT_TAXONOMY_ENDPOINT
        all_descendants = set()  # Using a set to avoid duplicates

        # We need to get all children at the requested level, AND at lower levels. That's what we're using the ranks array
        # for.
        rank_idx = TaxonManager.NCBI_RANKS.index(descendants_rank)
        descendants_ranks = TaxonManager.NCBI_RANKS[rank_idx:]

        # Split the target taxa into batches of 15
        for i in range(0, len(target_taxa), TaxonManager.TAXONOMY_ENDPOINT_BATCH_SIZE):
            batch = target_taxa[i:i+TaxonManager.TAXONOMY_ENDPOINT_BATCH_SIZE]

            # Prepare the request payload
            payload = {
                "input": batch,
                "descendants": True,
                "descendants_ranks": descendants_ranks
            }

            # Perform the HTTP POST request
            response = RequestManager.perform_post_request(url, payload)

            # Check if the request was successful
            if response.status_code == 200:
                data = response.json()

                # Extract descendants from each item in the response
                for item in data:
                    all_descendants.update(item.get("descendants", []))
            else:
                logging.error(f"Failed to retrieve taxonomy data for batch {batch}. Status code: {response.status_code}")

        # Convert the set of descendants back to a list and return it
        return list(all_descendants)


    @staticmethod
    def get_lineages_for_taxa(target_taxa: List[int]) -> Dict[int, List[int | None]]:
        url = TaxonManager.UNIPEPT_URL + TaxonManager.UNIPEPT_TAXONOMY_ENDPOINT

        # Remove duplicates from input and filter those already cached
        target_taxa = set(target_taxa)
        lineages = dict()

        # Prepare a list of taxa that are not yet in the cache
        taxa_to_request = [taxon for taxon in target_taxa if taxon not in TaxonManager.lineage_cache]

        # Fetch lineages from the API for taxa not in the cache
        for i in range(0, len(taxa_to_request), TaxonManager.TAXONOMY_ENDPOINT_BATCH_SIZE):
            batch = taxa_to_request[i:i+TaxonManager.TAXONOMY_ENDPOINT_BATCH_SIZE]

            payload = {
                "input": batch,
                "extra": True
            }

            # Perform the HTTP POST request
            response = RequestManager.perform_post_request(url, payload)

            if response.status_code == 200:
                data = response.json()
                for item in data:
                    lineage = [item.get(rank + "_id") for rank in TaxonManager.NCBI_RANKS]
                    taxon_id = item["taxon_id"]
                    # Cache the retrieved lineage
                    TaxonManager.lineage_cache[taxon_id] = lineage
            else:
                logging.error(f"Failed to retrieve taxonomy data for batch {batch}. Status code: {response.status_code}")

        # Collect results from the cache for all requested taxa
        for taxon in target_taxa:
            lineages[taxon] = TaxonManager.lineage_cache.get(taxon)

        return lineages


    @staticmethod
    def get_names_for_taxa(target_taxa: List[int]) -> Dict[int, str]:
        """
        Returns a mapping from taxon ID to taxon name for all taxa that have been provided to this function.

        Parameters
        ----------
        target_taxa: List[int],
            A list of taxon IDs for which all corresponding taxon names should be retrieved.
        """
        url = TaxonManager.UNIPEPT_URL + TaxonManager.UNIPEPT_TAXONOMY_ENDPOINT

        output = dict()

        # Split the target taxa into batches of 15
        for i in range(0, len(target_taxa), TaxonManager.TAXONOMY_ENDPOINT_BATCH_SIZE):
            batch = target_taxa[i:i+TaxonManager.TAXONOMY_ENDPOINT_BATCH_SIZE]

            # Prepare the request payload
            payload = {
                "input": batch
            }

            # Perform the HTTP POST request
            response = RequestManager.perform_post_request(url, payload)

            # Check if the request was successful
            if response.status_code == 200:
                data = response.json()

                # Extract descendants from each item in the response
                for item in data:
                    output[item.get("taxon_id")] = item.get("taxon_name")
            else:
                logging.error(f"Failed to retrieve taxonomy data for batch {batch}. Status code: {response.status_code}")

        return output