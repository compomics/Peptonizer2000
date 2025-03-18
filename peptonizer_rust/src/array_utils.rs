
pub fn normalize(mut array: Vec<f64>) -> Vec<f64> {
    let sum: f64 = array.iter().sum();
    for val in array.iter_mut() {
        *val /= sum;
    }
    array
}

pub fn log_normalize(mut array: Vec<f64>) -> Vec<f64> {
    let max_val = array.iter().cloned().fold(f64::NEG_INFINITY, f64::max); // Find the max value to prevent overflow
    let log_sum_exp = array.iter().map(|&x| (x - max_val).exp()).sum::<f64>().ln(); // Calculate logsumexp

    array.iter_mut()
        .map(|&mut x| (x - max_val - log_sum_exp).exp()) // Log-normalize and apply exp to each element
        .collect()
}

pub fn avoid_underflow(mut array: Vec<f64>) -> Vec<f64> {
    array.iter_mut().for_each(|x| if *x < 1e-30 { *x = 1e-30 });

    array
}