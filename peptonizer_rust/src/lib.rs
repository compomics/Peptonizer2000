use wasm_bindgen::prelude::*;
use serde_json::{Value, Result};
use utils::*;
use std::collections::{HashMap, HashSet};
use random::select_random_samples_with_weights;


extern crate wasm_bindgen;
extern crate serde_json;
extern crate serde;
extern crate web_sys;
extern crate wasm_bindgen_futures;
extern crate js_sys;

mod taxon_manager;
mod utils;
mod http_client;
mod random;

#[wasm_bindgen]
pub fn perform_taxa_weighing(
    unipept_responses: String,
    pep_scores: String,
    pep_psm_counts: String,
    max_taxa: usize,
    taxa_rank: String
) {
    log("Parsing Unipept responses from disk...");
    let unipept_responses: Vec<UnipeptJson> = serde_json::from_str(&unipept_responses).unwrap();

    let sequences: Vec<String> = unipept_responses.iter().map(|response| response.sequence.to_owned()).collect();

    let mut taxa: Vec<Vec<i32>> = unipept_responses.into_iter().map(|response| response.taxa).collect();
    
    log("Started mapping all taxon ids to the specified rank...");
    normalize_unipept_responses(&mut taxa, &taxa_rank);
    let chosen_idx: Vec<usize> = weighted_random_sample(&sequences, &taxa, 15000);

    log(&format!("Using {} sequences as input...", chosen_idx.len()));

    log("Normalizing peptides and converting to vectors...");

    let sequences: Vec<String> = chosen_idx.iter().map(|idx| sequences[*idx].to_owned()).collect();
    let taxa: Vec<Vec<i32>> = chosen_idx.iter().map(|idx| taxa[*idx].to_owned()).collect();

    // Parse scores from JSON string to hashmap, only keep the randomly selected samples.
    let pep_scores_map: HashMap<String, f32> = serde_json::from_str(&pep_scores).unwrap();
    let mut pep_scores: Vec<f32> = vec![0.0; sequences.len()];
    for i in 0..sequences.len() {
        pep_scores[i] = pep_scores_map[&sequences[i]];
    }

    // parse counts from JSON string to hashmap, only keep the randomly selected samples.
    let pep_psm_counts_map: HashMap<String, i32> = serde_json::from_str(&pep_psm_counts).unwrap();
    let mut pep_psm_counts: Vec<i32> = vec![0; sequences.len()];
    for i in 0..sequences.len() {
        pep_psm_counts[i] = pep_psm_counts_map[&sequences[i]];
    }

    /* Score the degeneracy of a taxa, i.e.,
       how conserved a peptide sequence is between taxa.
       map all taxids in the list in the taxa column back to their taxid at species level (or the rank specified by the user)
       Right now, HigherTaxa is simply a copy of taxa. This step still needs to be optimized.
       Move taxa to highertaxa because taxa is not used anymore.
    */
    let higher_taxa: Vec<Vec<i32>> = taxa; 

    // Divide the number of PSMs of a peptide by the number of taxa the peptide is associated with, exponentiated by 3
    log("Started dividing the number of PSMS of a peptide by the number the peptide is associated with...");
    let weights: Vec<f32> = pep_psm_counts.iter()
                                          .zip(higher_taxa.iter().map(|taxa| taxa.len().pow(3)))
                                          .map(|(&count, len_cube)| count as f32 / len_cube as f32)
                                          .collect();

    let unique_psm_taxa: HashSet<i32> = higher_taxa.iter()
                                                   .filter(|tax| tax.len() == 1)
                                                   .map(|tax| tax[0])
                                                   .collect();

    // Sum up the weights of a taxon and sort by weight
    log("Started summing the weights of a taxon and sorting them by weight...");
    let log_weight: Vec<f32> = weights.iter().map(|w| (w + 1.0).log10()).collect();

    let mut tax_id_weights: HashMap<i32, f32> = HashMap::new();
    for (ids, weight) in higher_taxa.into_iter().zip(log_weight.clone().into_iter()) {
        for id in ids {
            *tax_id_weights.entry(id).or_insert(0.0) += weight;
        }
    }
    let mut sorted_tax_id_weights: Vec<(i32, f32)> = tax_id_weights.into_iter().collect();
    sorted_tax_id_weights.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
    let (tax_ids, tax_id_weights): (Vec<i32>, Vec<f32>) = sorted_tax_id_weights.into_iter().unzip();

    //  Since large proteomes tend to have more detectable peptides,
    // we adjust the weight by dividing by the size of the proteome i.e.,
    // the number of proteins that are associated with a taxon
    let scaled_weight = log_weight;

    // Retrieves the specified taxonomic rank taxid in the lineage of each of the species-level taxids returned by
    // Unipept for both the UnipeptFrame and the TaxIdWeightFrame
    let higher_unique_psm_taxids = unique_psm_taxa;

    // Group the duplicate entries of higher up taxa and sum their weights
    let higher_taxid_weights = scaled_weight;

    let tax_ids: Vec<i32> = tax_ids.into_iter().filter(|id| *id != 1869227).collect(); // TODO: is this correct?
    let higher_taxid_unique = tax_ids.iter().map(|id| higher_unique_psm_taxids.contains(&id));

    // TODO: check if duplicate removal is necessary, example datasets do not contain doubles.

    if tax_ids.len() < 50 {
        // Return csvs
        log("max taxa < 50");
    } else {
        let mut taxa_to_include: HashSet<i32> = tax_ids.iter().take(max_taxa).cloned().collect();
        taxa_to_include.extend(higher_unique_psm_taxids);
        // filter score vectors and return csvs
    }
}

fn normalize_unipept_responses(taxa: &mut Vec<Vec<i32>>, taxa_rank: &str) {
    
    let mut lineage_cache: HashMap<i32, Vec<Option<i32>>> = HashMap::new();

    // Map all taxa onto the rank specified by the user
    for i in 0..taxa.len() {
        taxa[i] = taxon_manager::TaxonManager::get_unique_lineage_at_specified_rank(&taxa[i], taxa_rank, &mut lineage_cache);
    }
}

fn weighted_random_sample(sequences: &Vec<String>, taxa: &Vec<Vec<i32>>, n: usize) -> Vec<usize> {
    
    // Calculate normalized weights based on the length of the taxa array
    let weights: Vec<f64> = taxa.iter().map(|taxon| if taxon.len() == 0 { 0.0 } else { 1.0 / taxon.len() as f64 }).collect();
    let total_weight: f64 = weights.iter().sum();
    let normalized_weights: Vec<f64> = weights.iter().map(|w| w / total_weight as f64).collect();

    let samples: Vec<usize> = select_random_samples_with_weights(sequences, normalized_weights, n);

    samples
}
