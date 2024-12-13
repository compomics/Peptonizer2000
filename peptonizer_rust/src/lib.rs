use wasm_bindgen::prelude::*;
use serde::Deserialize;
use serde_json::{Value, Result};
use utils::*;
use std::collections::HashMap;

extern crate wasm_bindgen;
extern crate serde_json;
extern crate serde;
extern crate web_sys;
extern crate wasm_bindgen_futures;

mod taxon_manager;
mod utils;
mod http_client;

#[wasm_bindgen]
pub fn perform_taxa_weighing(
    unipept_responses: String,
    pep_scores: String,
    pep_psm_counts: String,
    max_taxa: i32,
    taxa_rank: String
) {
    log("Parsing Unipept responses, scores and counts from disk...");
    let mut unipept_responses: Vec<UnipeptJson> = serde_json::from_str(&unipept_responses).unwrap();

    let pep_scores_map: HashMap<String, f32> = serde_json::from_str(&pep_scores).unwrap();
    let mut pep_scores: Vec<f32> = vec![0.0; unipept_responses.len()];
    for i in 0..unipept_responses.len() {
        pep_scores[i] = pep_scores_map[&unipept_responses[i].sequence];
    }

    let pep_psm_counts_map: HashMap<String, i32> = serde_json::from_str(&pep_psm_counts).unwrap();
    let mut pep_psm_counts: Vec<i32> = vec![0; unipept_responses.len()];
    for i in 0..unipept_responses.len() {
        pep_psm_counts[i] = pep_psm_counts_map[&unipept_responses[i].sequence];
    }

    log("Started mapping all taxon ids to the specified rank...");
    normalize_unipept_responses(&mut unipept_responses, &taxa_rank);
}

#[derive(Debug, Deserialize)]
struct UnipeptJson {
    sequence: String,
    taxa: Vec<i32>,
}

fn normalize_unipept_responses(unipept_responses: &mut Vec<UnipeptJson>, taxa_rank: &str) {
    
    let mut lineage_cache: HashMap<i32, Vec<Option<i32>>> = HashMap::new();

    // Map all taxa onto the rank specified by the user
    for i in 0..unipept_responses.len() {
        unipept_responses[i].taxa = get_unique_lineage_at_specified_rank(&unipept_responses[i].taxa, taxa_rank, &mut lineage_cache);
    }
}

fn get_unique_lineage_at_specified_rank(tax_ids: &Vec<i32>, taxa_rank: &str, lineage_cache: &mut HashMap<i32, Vec<Option<i32>>>) -> Vec<i32> {
    
    let lineages = taxon_manager::TaxonManager::get_lineages_for_taxa(tax_ids, lineage_cache);
    
    return vec!(1);
}
