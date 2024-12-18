#[cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;

#[cfg(not(target_arch = "wasm32"))]
use pyo3::prelude::*;

use crate::weight_taxa::perform_taxa_weighing;

extern crate serde_json;
extern crate serde;

mod taxon_manager;
mod utils;
mod http_client;
mod random;
mod weight_taxa;

#[cfg(target_arch = "wasm32")]
mod wasm_deps {
    extern crate wasm_bindgen;
    extern crate web_sys;
    extern crate wasm_bindgen_futures;
    extern crate js_sys;
}

#[cfg(not(target_arch = "wasm32"))]
mod py_deps {

}

#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
pub fn perform_taxa_weighing_wasm(
    unipept_responses: String,
    pep_scores: String,
    pep_psm_counts: String,
    max_taxa: usize,
    taxa_rank: String
) -> Box<[JsValue]> {
    
    let (sequence_csv, taxa_weights_csv): (String, String) = perform_taxa_weighing(unipept_responses, pep_scores, pep_psm_counts, max_taxa, taxa_rank);

    Box::new([JsValue::from(sequence_csv), JsValue::from(taxa_weights_csv)])
}


#[cfg(not(target_arch = "wasm32"))]
#[pyfunction]
fn perform_taxa_weighing_py(
    unipept_responses: String,
    pep_scores: String,
    pep_psm_counts: String,
    max_taxa: usize,
    taxa_rank: String
) -> (String, String) {
    perform_taxa_weighing(unipept_responses, pep_scores, pep_psm_counts, max_taxa, taxa_rank)
}

#[cfg(not(target_arch = "wasm32"))]
#[pymodule]
fn peptonizer_rust(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(perform_taxa_weighing_py, m)?)?;
    Ok(())
}