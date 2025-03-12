use crate::factor_graph::CTFactorGraph;
use crate::node::{Node, NodeType};
use std::collections::{HashMap, HashSet};
use std::mem;
use crate::utils::log;

#[derive(Debug, Clone)]
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
    graph: CTFactorGraph,
    max_val: Option<(i32, i32)>,
    priorities: HashMap<i32, f64>,
    // Keeps track of residuals for duos of directed edges, indexed by [start node][id of end node in start neighbours][id of neighbour in end node]
    total_residuals: Vec<Vec<Vec<f64>>>, 
    // Maps a node ID onto its current belief value
    current_beliefs: Vec<NodeBelief>, 
    // incoming messages for each node [end node][neighbour id]
    msg_in: Vec<Vec<Vec<f64>>>,
    msg_in_new: Vec<Vec<Vec<f64>>>,
    msg_in_log: Vec<Vec<Vec<f64>>>
}

impl Messages {

    pub fn new(ct_graph_in: CTFactorGraph) -> Messages {
        let max_val: Option<(i32, i32)> = None;
        
        // TODO
        let priorities = HashMap::new();

        let mut total_residuals: Vec<Vec<Vec<f64>>> = Vec::with_capacity(ct_graph_in.node_count());
        for node in ct_graph_in.get_nodes() {
            
            let mut total_residual_node: Vec<Vec<f64>> = Vec::with_capacity(node.neighbors_count());
            for neighbor in ct_graph_in.get_neighbors(node) {
                total_residual_node.push(vec![0.0; ct_graph_in.get_node(neighbor).neighbors_count()]);
            }

            total_residuals.push(total_residual_node);

        }

        let mut current_beliefs: Vec<NodeBelief> = Vec::with_capacity(ct_graph_in.node_count());
        for node in ct_graph_in.get_nodes() {
            current_beliefs.push(get_initial_belief(node));
        }

        let mut msg_in = Vec::with_capacity(ct_graph_in.node_count());
        let mut msg_in_new = Vec::with_capacity(ct_graph_in.node_count());
        for node in ct_graph_in.get_nodes() {
            
            let mut msg_in_node = Vec::with_capacity(node.neighbors_count());
            let mut msg_in_new_node = Vec::with_capacity(node.neighbors_count());
            for edge_id in node.get_incident_edges() {
                match ct_graph_in.get_edge(*edge_id).get_message_length() {
                    Some(message_length) => {
                        msg_in_node.push(vec![1.0; message_length as usize]);
                        msg_in_new_node.push(vec![1.0; message_length as usize]);
                    },
                    None => {
                        msg_in_node.push(vec![0.5, 0.5]);
                        msg_in_new_node.push(vec![0.0, 0.0]);
                    }
                }
            }

            msg_in.push(msg_in_node);
            msg_in_new.push(msg_in_new_node);
        }

        let msg_in_log = msg_in.clone();

        Messages { graph: ct_graph_in, max_val, priorities, total_residuals, current_beliefs, msg_in, msg_in_new, msg_in_log }

    }

    
}