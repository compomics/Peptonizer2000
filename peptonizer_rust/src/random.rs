use js_sys::Math;
use std::collections::HashSet;

#[cfg(target_arch = "wasm32")]
pub fn select_random_samples_with_weights(
    weights: Vec<f64>,
    n: usize,
) -> Vec<usize> {

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

    let samples: Vec<usize> = samples.iter().cloned().collect();

    samples
}