use crate::factor_graph::CTFactorGraph;
use crate::node::{Node, NodeType};
use std::collections::HashMap;

#[derive(Debug)]
pub enum NodeBelief {
    PeptideBelief(f64, f64),
    FactorBelief(Vec<[f64;2]>),
    TaxonBelief(f64, f64),
    ConvolutionTreeBelief
}

pub fn get_initial_belief(node: &Node) -> NodeBelief {
    match node.get_subtype() {
        NodeType::PeptideNode { initial_belief_0, initial_belief_1 } => NodeBelief::PeptideBelief(*initial_belief_0, *initial_belief_1),
        NodeType::FactorNode { initial_belief, .. } => NodeBelief::FactorBelief(initial_belief.array.clone()),
        NodeType::TaxonNode { initial_belief_0, initial_belief_1 } => NodeBelief::TaxonBelief(*initial_belief_0, *initial_belief_1),
        NodeType::ConvolutionTreeNode { .. } => NodeBelief::ConvolutionTreeBelief
    }
}

pub struct Messages {
    max_val: Option<(i32, i32)>,
    priorities: HashMap<i32, f64>,
    total_residuals: Vec<Vec<f64>>, // Keeps track of residuals for duos of edges, indexed by [edge index][neighbour index of edge's end node]
    current_beliefs: Vec<NodeBelief>, // Maps a node ID onto its current belief value
    msg_in: Vec<Vec<Vec<f64>>>,
    msg_in_new: Vec<Vec<Vec<f64>>>,
    msg_in_log: Vec<Vec<Vec<f64>>>
}

impl Messages {

    pub fn new(ct_graph_in: CTFactorGraph) {
        let mut total_residuals: Vec<Vec<f64>> = Vec::with_capacity(ct_graph_in.edge_count());
        for edge in ct_graph_in.get_edges() {
            let length = ct_graph_in.get_node(edge.get_node2_id()).neighbors_count();
            total_residuals.push(vec![0.0; length]);
        }

        let mut current_beliefs: Vec<NodeBelief> = Vec::with_capacity(ct_graph_in.node_count());
        for node in ct_graph_in.get_nodes() {
            current_beliefs.push(get_initial_belief(node));
        }

        
    }
}