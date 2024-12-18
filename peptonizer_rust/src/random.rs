#[cfg(target_arch = "wasm32")]
use js_sys::Math;

use std::collections::HashSet;

#[cfg(target_arch = "wasm32")]
pub fn select_random_samples_with_weights(
    weights: Vec<f64>,
    n: usize,
) -> HashSet<usize> {

    let cumulative_weights: Vec<f64> = weights
        .iter()
        .scan(0.0, |acc, &weight| {
            *acc += weight;
            Some(*acc)
        })
        .collect();

    let mut samples: HashSet<usize> = HashSet::with_capacity(n);

    while samples.len() < n {
        let r = Math::random() as f64; 
        if r >= cumulative_weights[cumulative_weights.len()-1] {
            continue;
        }

        let chosen_idx: usize = cumulative_weights.binary_search_by(|&w| w.partial_cmp(&r).unwrap()).unwrap_or_else(|x| x);

        samples.insert(chosen_idx);
    }

    samples
}

#[cfg(not(target_arch = "wasm32"))]
pub fn select_random_samples_with_weights(
    weights: Vec<f64>,
    n: usize,
) -> HashSet<usize> {
    use rand::prelude::*;
    use rand::distributions::WeightedIndex;

    // Create a WeightedIndex for sampling
    let dist = WeightedIndex::new(&weights).unwrap();

    let mut rng = thread_rng();
    let mut selected_indices = HashSet::new();

    while selected_indices.len() < n {
        let idx = dist.sample(&mut rng);
        selected_indices.insert(idx);
    }

    selected_indices
}
