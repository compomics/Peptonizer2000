use crate::factor_graph::CTFactorGraph;

/// Runs the belief propagation algorithm on a graph that's represented by the string in graphml_content with the
/// tuning parameters further specified to this function. This function returns a string that contains the result of
/// the belief propagation algorithm, represented as a CSV (and can thus directly be written to a CSV-file, if desired).
pub fn run_belief_propagation(
    graph: String,
    alpha: f32,
    beta: f32,
    regularized: bool,
    prior: f32,
    max_iter: i32,
    tol: f32
) -> String {
    
    let graph = CTFactorGraph::from_graphml(&graph).unwrap();

    format!("{}", graph.to_string())
}