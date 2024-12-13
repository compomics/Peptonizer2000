use js_sys::Math;
use crate::utils::*;

#[cfg(target_arch = "wasm32")]
pub fn select_random_responses_with_weights(
    responses: Vec<UnipeptJson>,
    weights: Vec<f32>,
    n: usize,
) -> Vec<UnipeptJson> {

    let cumulative_weights: Vec<f32> = weights
        .iter()
        .scan(0.0, |acc, &weight| {
            *acc += weight;
            Some(*acc)
        })
        .collect();

    let mut samples: Vec<UnipeptJson> = Vec::with_capacity(n);

    log(&format!("start {}", cumulative_weights[cumulative_weights.len()-1]));
    while samples.len() < n {
        let r = Math::random() as f32; 
        if r >= cumulative_weights[cumulative_weights.len()-1] {
            continue;
        }

        let mut chosen_idx: usize = cumulative_weights.binary_search_by(|&w| w.partial_cmp(&r).unwrap()).unwrap_or_else(|x| x);
        samples.push(responses[chosen_idx].clone());
    }

    samples
}