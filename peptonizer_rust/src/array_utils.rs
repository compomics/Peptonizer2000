
pub fn normalize(array: &mut Vec<f64>) {
    let sum: f64 = array.iter().sum();
    for val in array.iter_mut() {
        *val /= sum;
    }
}

pub fn normalize_2d(array: &mut Vec<[f64; 2]>) {
    let sum: f64 = array.iter().map(|x| x[0] + x[1]).sum();
    for val in array.iter_mut() {
        val[0] /= sum;
        val[1] /= sum;
    }
}

pub fn log_normalize(array: &mut Vec<f64>) {
    let max_val = array.iter().cloned().fold(f64::NEG_INFINITY, f64::max); // Find the max value to prevent overflow
    let log_sum_exp = array.iter().map(|&x| (x - max_val).exp()).sum::<f64>().ln(); // Calculate logsumexp

    array.iter_mut()
        .for_each(|x| *x = (*x - max_val - log_sum_exp).exp()); // Log-normalize and apply exp to each element
}

pub fn log_normalize_2d(array: &mut Vec<[f64;2]>) {
    let max_val = array.iter().cloned().fold(f64::NEG_INFINITY, |acc, [a, b]| f64::max(f64::max(a, b), acc)); // Find the max value to prevent overflow
    let log_sum_exp = array.iter().map(|&[a, b]| (a - max_val).exp() + (b - max_val).exp()).sum::<f64>().ln(); // Calculate logsumexp

    array.iter_mut()
        .for_each(|x| {
            x[0] = (x[0] - max_val - log_sum_exp).exp();
            x[1] = (x[1] - max_val - log_sum_exp).exp();
        }); // Log-normalize and apply exp to each element
}

pub fn avoid_underflow(array: &mut Vec<f64>) {
    array.iter_mut().for_each(|x| if *x < 1e-30 { *x = 1e-30 });
}