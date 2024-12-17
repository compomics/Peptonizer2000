use std::collections::{ HashMap, HashSet };
use serde_json::{ Value };
use crate::http_client::*;

pub struct TaxonManager;

impl TaxonManager{
    const UNIPEPT_URL: &str = "http://api.unipept.ugent.be";
    const UNIPEPT_TAXONOMY_ENDPOINT: &str = "/api/v2/taxonomy.json";

    const TAXONOMY_ENDPOINT_BATCH_SIZE: usize = 100;

    const NCBI_RANKS: &[&str] = &[
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
    ];

    fn parse_response_json_string(http_response: &str) -> Vec<HashMap<String, Option<i32>>> {
        let http_response_map: Vec<HashMap<String, Option<i32>>> = serde_json::from_str::<Vec<HashMap<String, Value>>>(http_response)
                .unwrap()
                .into_iter()
                .map(|mut obj: HashMap<String, Value>| {
                    // Remove the key-value pair where the value is a string
                    obj.retain(|_, v| v.is_null() || v.is_number());

                    // Convert the remaining keys and values to `HashMap<String, Option<i32>>`
                    obj.into_iter()
                        .map(|(key, value)| {
                            let value = if value.is_number() {
                                Some(value.as_i64().unwrap() as i32)
                            } else {
                                None
                            };
                            (key, value)
                        })
                        .collect()
                })
                .collect();
        
        http_response_map
    }

    pub fn get_unique_lineage_at_specified_rank(target_taxa: &Vec<i32>, taxa_rank: &str, lineage_cache: &mut HashMap<i32, Vec<Option<i32>>>) -> Vec<i32> {

        let url: String = format!("{}{}", TaxonManager::UNIPEPT_URL, TaxonManager::UNIPEPT_TAXONOMY_ENDPOINT);

        // Remove duplicates from input
        let target_taxa: HashSet<i32> = target_taxa.iter().cloned().collect();

        // Prepare a list of taxa that are not yet in the cache
        let taxa_to_request: Vec<i32> = target_taxa.iter().filter(| tax_id | ! lineage_cache.contains_key(tax_id)).cloned().collect();

        // Fetch lineages from the API for taxa not in the cache
        for i in (0..taxa_to_request.len()).step_by(TaxonManager::TAXONOMY_ENDPOINT_BATCH_SIZE) {

            let batch_size: usize = std::cmp::min(TaxonManager::TAXONOMY_ENDPOINT_BATCH_SIZE, taxa_to_request.len() - i);
            let batch: Vec<i32> = taxa_to_request[i..(i + batch_size)].to_vec();

            // Perform the HTTP POST request
            let http_client: &dyn HttpClient = &create_http_client();
            let http_response:  String = 
                http_client.perform_post_request(url.clone(), batch)
                .map_err(|e| format!("Failed to retrieve taxonomy data for batch {}. Error message: {}", (i / TaxonManager::TAXONOMY_ENDPOINT_BATCH_SIZE), e))
                .unwrap();
            let http_response = Self::parse_response_json_string(&http_response);

            for lineage_json in &http_response {
                let lineage: Vec<Option<i32>> = TaxonManager::NCBI_RANKS.iter()
                        .filter_map(|key| lineage_json.get(&format!("{}_id", key)).cloned())
                        .collect();
                let taxon_id: i32 = lineage_json.get("taxon_id").unwrap().unwrap();
                lineage_cache.insert(taxon_id, lineage);
            }
            
        }

        let rank_idx = TaxonManager::NCBI_RANKS.iter().position(|&ncbi_rank| ncbi_rank == taxa_rank).unwrap();
        let lineage: HashSet<i32> = target_taxa.iter()
                                                .filter_map(|taxon| lineage_cache[&taxon][rank_idx].clone())
                                                .collect();
        let lineage: Vec<i32> = lineage.into_iter().collect();

        lineage
    }

}
        