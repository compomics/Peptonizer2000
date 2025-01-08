use crate::factor_graph::CTFactorGraph;
use std::collections::HashMap;


/// Performs bayesian inference through loopy belief propagation, returns dictionary {variable:posterior_probability}
///
/// # Arguments
/// * ct_factor_graphs: list, contains FactorGraph objects on which inference can be performed
/// * max_iterations: int, max number of iterations in case of non-convergence
/// * tolerance: float, error tolerance between messages for convergence criterion
fn calibrate_all_subgraphs(
    ct_factor_graphs: Vec<CTFactorGraph>,
    max_iterations: i32,
    tolerance: f64
) {
    let mut results: HashMap<String, Vec<[f64; 2]>> = HashMap::new();
    let mut node_categories: HashMap<String, String> = HashMap::new();

    for subgraph in ct_factor_graphs {
        if subgraph.node_count() > 2 {

            subgraph.add_node_categories(&mut node_categories);
            
        }
    }
}

/// Runs the belief propagation algorithm on a graph that's represented by the string in graphml_content with the
/// tuning parameters further specified to this function. This function returns a string that contains the result of
/// the belief propagation algorithm, represented as a CSV (and can thus directly be written to a CSV-file, if desired).
pub fn run_belief_propagation(
    graph: String,
    alpha: f64,
    beta: f64,
    regularized: bool,
    prior: f64,
    max_iter: i32,
    tol: f64
) -> String {
    
    let mut ct_factor_graph = CTFactorGraph::from_graphml(&graph).unwrap();
    ct_factor_graph.fill_in_factors(alpha, beta, regularized);
    ct_factor_graph.fill_in_priors(prior);
    ct_factor_graph.add_ct_nodes();

    let ct_factor_graphs = ct_factor_graph.connected_components();

    calibrate_all_subgraphs(
        ct_factor_graphs,
        max_iter,
        tol
    );

    format!("Hello")
}