
extern crate serde_json;
extern crate serde;

mod taxon_manager;
mod utils;
mod http_client;
mod random;
mod weight_taxa;
mod zero_lookahead_belief_propagation;
mod node;
mod factor_graph;
mod messages;

#[cfg(target_arch = "wasm32")]
pub use wasm::*;

#[cfg(not(target_arch = "wasm32"))]
pub use pyo3::*;

#[cfg(target_arch = "wasm32")]
mod wasm {
    use wasm_bindgen::prelude::*;
    use crate::weight_taxa::perform_taxa_weighing;
    use crate::zero_lookahead_belief_propagation::run_belief_propagation;
    use crate::utils::log;

    extern crate wasm_bindgen;
    extern crate web_sys;
    extern crate wasm_bindgen_futures;
    extern crate js_sys;

    #[wasm_bindgen]
    pub fn perform_taxa_weighing_wasm(
        unipept_responses: String,
        pep_scores: String,
        pep_psm_counts: String,
        max_taxa: usize,
        taxa_rank: String
    ) -> Box<[JsValue]> {
        log("1");
        let (sequence_csv, taxa_weights_csv): (String, String) = perform_taxa_weighing(unipept_responses, pep_scores, pep_psm_counts, max_taxa, taxa_rank);
        log("2");
        Box::new([JsValue::from(sequence_csv), JsValue::from(taxa_weights_csv)])
    }

    #[wasm_bindgen]
    pub fn run_belief_propagation_wasm(
        graph: String,
        alpha: f64,
        beta: f64,
        regularized: bool,
        prior: f64,
        max_iter: Option<i32>,
        tol: Option<f64>
    ) -> String {
        let max_iter: i32 = max_iter.unwrap_or(10000);
        let tol: f64 = tol.unwrap_or(0.006);
        
        run_belief_propagation(graph, alpha, beta, regularized, prior, max_iter, tol)
    }

}

#[cfg(not(target_arch = "wasm32"))]
mod pyo3 {
    use pyo3::prelude::*;
    use crate::weight_taxa::perform_taxa_weighing;

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
}
