use crate::factor_graph::CTFactorGraph;
use crate::node::{Node, NodeType};
use std::collections::{HashMap, HashSet};
use std::mem;
use crate::utils::log;
use crate::array_utils::*;
use crate::convolution_tree::ConvolutionTree;

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

impl NodeBelief {
    pub fn values(&self) -> Vec<f64> {
        match self {
            NodeBelief::PeptideBelief(a, b) => vec![*a, *b],
            NodeBelief::FactorBelief(vec) => vec.iter().flat_map(|arr| arr.to_vec()).collect(),
            NodeBelief::TaxonBelief(a, b) => vec![*a, *b],
            NodeBelief::ConvolutionTreeBelief => vec![],
        }
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

    pub fn zero_lookahead_bp(&mut self, max_loops: i32, tolerance: f64) -> Vec<NodeBelief> {

        let mut max_residual: f64 = f64::MAX;

        // first, do 5 loops where I update all messages
        for _ in 0..5 {
            self.compute_update(false);
            let temp = mem::replace(&mut self.msg_in_log, mem::take(&mut self.msg_in));
            self.msg_in = mem::take(&mut self.msg_in_new);
            self.msg_in_new = temp;
        }

        self.current_beliefs.clone()
    }

    fn compute_update(&mut self, local_loops: bool) {
        let mut checked_cts: HashSet<i32> = HashSet::new();
        
        if self.max_val.is_some() && local_loops {

            let (_, start_id) = self.max_val.unwrap();
            let start_node = self.graph.get_node(start_id);
            for (end_in_start_id, end_id) in self.graph.get_neighbors(start_node).iter().enumerate() {
                self.single_message_update(start_id, *end_id, end_in_start_id as i32, &mut checked_cts);
            }

        } else {

            for id in 0..self.graph.node_count() {
                let id = id as i32;
                let start_node = self.graph.get_node(id);
                for (end_in_start_id, end_id) in self.graph.get_neighbors(start_node).iter().enumerate() {
                    self.single_message_update(id, *end_id, end_in_start_id as i32, &mut checked_cts);
                }
            }

        }
    }

    fn single_message_update(&mut self, start_id: i32, end_id: i32, end_in_start_id: i32, checked_cts: &mut HashSet<i32>) {

        let start_node: &Node = self.graph.get_node(start_id);

        match start_node.get_subtype() {
            NodeType::TaxonNode { .. } | NodeType::PeptideNode { .. } => {
                let new_message = self.compute_out_message_variable(start_id, end_id, end_in_start_id);
                let end_node: &Node = self.graph.get_node(end_id);
                let start_in_end_id = self.graph.get_neighbors(end_node).iter().position(|neighbor_id| *neighbor_id == start_id).unwrap();
                self.msg_in_new[end_id as usize][start_in_end_id] = new_message;
            },
            NodeType::FactorNode { .. } => {
                let new_message = self.compute_out_message_factor(start_id, end_id, end_in_start_id);
                let end_node: &Node = self.graph.get_node(end_id);
                let start_in_end_id = self.graph.get_neighbors(end_node).iter().position(|neighbor_id| *neighbor_id == start_id).unwrap();
                self.msg_in_new[end_id as usize][start_in_end_id] = new_message;
            },
            NodeType::ConvolutionTreeNode { .. } => 
                if ! checked_cts.contains(&start_id) {
                    self.compute_out_messages_ct_tree(start_id);
                    checked_cts.insert(start_id);
                }
        };
    }

    fn compute_out_message_variable(&mut self, start_id: i32, end_id: i32, end_in_start_id: i32) -> Vec<f64> {
        
        let start_node = self.graph.get_node(start_id);
        let end_node = self.graph.get_node(end_id);
        let mut incoming_messages_end: Vec<Vec<f64>> = self.msg_in[start_id as usize].clone();
        incoming_messages_end.remove(end_in_start_id as usize);
        let node_belief: Vec<f64> = self.current_beliefs[start_id as usize].values();

        if incoming_messages_end.len() == 0 {
            // TODO: check in Python code the idea of the any statement (ask Tanja)
            let start_in_end_id = self.graph.get_neighbors(end_node).iter().position(|neighbor_id| *neighbor_id == start_id).unwrap();
            return self.msg_in[end_id as usize][start_in_end_id as usize].clone();
        }

        // need for logs to prevent underflow in very large multiplications
        let incoming_messages_end: Vec<Vec<f64>> = incoming_messages_end
        .iter()
        .map(|msg| msg.iter().map(|&x| x.ln()).collect())
        .collect();

        // TODO: remove if this doesn't print
        for i in 0..incoming_messages_end.len() {
            if incoming_messages_end[i].len() != 2 {
                log("WRONG LENGTH");
            }
        }

        // Sum of incoming messages
        let sum_logs: Vec<f64> = incoming_messages_end.iter().fold(vec![0.0;2], |mut acc,  row| {acc[0] += row[0]; acc[1] += row[1]; acc});

        // Log transform node_belief
        let node_belief_log: Vec<f64> = node_belief.iter().map(|&x| x.ln()).collect();

        // Compute final log-normalized message
        let mut out_message_log: Vec<f64> = log_normalize(
            node_belief_log.iter().zip(sum_logs.iter()).map(|(&a, &b)| a + b).collect()
        );

        // Prevent underflow: Replace zeros with 1e-30
        out_message_log.iter_mut().for_each(|x| if *x == 0.0 { *x = 1e-30 });

        out_message_log
    }

    fn compute_out_message_factor(&mut self, start_id: i32, end_id: i32, end_in_start_id: i32) -> Vec<f64> {
        let start_node = self.graph.get_node(start_id);
        let end_node = self.graph.get_node(end_id);
        let mut incoming_messages_end: Vec<Vec<f64>> = self.msg_in[start_id as usize].clone();
        incoming_messages_end.remove(end_in_start_id as usize);
        let node_belief: Vec<f64> = self.current_beliefs[start_id as usize].values();

        match end_node.get_subtype() {
            NodeType::ConvolutionTreeNode { .. } => {
                // handles empty & messages with only one value
                incoming_messages_end.push(vec![1.0, 1.0]);
                
                // Product of incoming messages
                let prod: Vec<f64> = incoming_messages_end.iter().fold(vec![1.0;2], |mut acc,  row| {acc[0] *= row[0]; acc[1] *= row[1]; acc});
                
                // Compute final normalized message
                let mut out_message: Vec<f64> = normalize(
                    node_belief.iter().zip(prod.iter()).map(|(&a, &b)| a * b).collect(),
                );

                return out_message;
                // TODO: check shapes of previous variables (compare with python)
            },
            _ => {
                if incoming_messages_end[0].len() > 2 {
                    // TODO check shapes
                    let incoming_messages_log: Vec<f64> = incoming_messages_end[0].iter().map(|&x| x.ln()).collect();
                    let node_belief_log: Vec<f64> = node_belief.iter().map(|&x| x.ln()).collect();

                    let mut out_message_log: Vec<f64> = log_normalize(
                        node_belief_log.iter().zip(incoming_messages_log.iter()).map(|(&a, &b)| a + b).collect()
                    );

                    // Prevent underflow: Replace zeros with 1e-30
                    out_message_log.iter_mut().for_each(|x| if *x == 0.0 { *x = 1e-30 });

                    return out_message_log;

                } else {

                    incoming_messages_end.push(vec![1.0, 1.0]);

                    let prod: Vec<f64> = incoming_messages_end.iter().fold(vec![1.0;2], |mut acc,  row| {acc[0] *= row[0]; acc[1] *= row[1]; acc});

                    // Compute final normalized message
                    let mut out_message: Vec<f64> = normalize(
                        node_belief.iter().zip(prod.iter()).map(|(&a, &b)| a * b).collect(),
                    );

                    return out_message;
                }
            }
        }
    }

    fn compute_out_messages_ct_tree(&mut self, start_id: i32) {
        let start_node = self.graph.get_node(start_id);

        let mut prot_prob_list: Vec<Vec<f64>> = Vec::new();
        let mut old_prot_prob_list: Vec<Vec<f64>> = Vec::new();

        let mut shared_likelihoods: Vec<f64> = Vec::new();
        let mut old_shared_likelihoods: Vec<f64> = Vec::new();

        let mut peptides: Vec<i32> = Vec::new();
        let mut prot_list: Vec<i32> = Vec::new();

        let mut last_neighbor_id: Option<usize> = None;
        for (neighbor_id_in_start, &neighbor_id) in self.graph.get_neighbors(start_node).iter().enumerate() {
            let neighbor = self.graph.get_node(neighbor_id);
            match neighbor.get_subtype() {
                
                NodeType::FactorNode { .. } => {
                    prot_prob_list.push(self.msg_in[start_id as usize][neighbor_id_in_start].clone());
                    old_prot_prob_list.push(self.msg_in[start_id as usize][neighbor_id_in_start].clone());
                    prot_list.push(neighbor_id);
                },
                _ => {
                    peptides.push(neighbor_id);
                    
                    // Prevent underflow: Replace zeros with 1e-30
                    shared_likelihoods.iter_mut().for_each(|x| if *x == 0.0 { *x = 1e-30 });

                    for (a, b) in shared_likelihoods.iter_mut().zip(self.msg_in_log[start_id as usize][neighbor_id_in_start].iter()) {
                        *a *= b;
                    }

                    last_neighbor_id = Some(neighbor_id_in_start);
                }
            }
        }

        if let Some(neighbor_id) = last_neighbor_id {
            // Prevent underflow: Replace zeros with 1e-30
            shared_likelihoods = avoid_underflow(shared_likelihoods);
            old_shared_likelihoods = shared_likelihoods.iter()
                .zip(self.msg_in_log[start_id as usize][neighbor_id].iter())
                .map(|(a, b)| a * b).collect();
        }

        // TODO: revision this
        if old_shared_likelihoods != shared_likelihoods && prot_prob_list != old_prot_prob_list {
            let convolution_tree = ConvolutionTree::new(shared_likelihoods, prot_prob_list);

            for (protein_id, protein) in prot_list.iter().enumerate() {
                let node_neighbor_index = self.graph.get_neighbors_from_id(*protein).iter().position(|&x| x == start_id).unwrap();
                self.msg_in_new[*protein as usize][node_neighbor_index] = avoid_underflow(convolution_tree.message_to_variable(protein_id));
            }

            for pep in peptides {
                let node_neighbor_index = self.graph.get_neighbors_from_id(pep).iter().position(|&x| x == start_id).unwrap();
                self.msg_in_new[pep as usize][node_neighbor_index] = avoid_underflow(convolution_tree.message_to_shared_likelihood());
            }

        } else {
            
            for protein in prot_list {
                let node_neighbor_index = self.graph.get_neighbors_from_id(protein).iter().position(|&x| x == start_id).unwrap();
                self.msg_in_new[protein as usize][node_neighbor_index] = self.msg_in[protein as usize][node_neighbor_index].clone();
            }

            for pep in peptides {
                let node_neighbor_index = self.graph.get_neighbors_from_id(pep).iter().position(|&x| x == start_id).unwrap();
                self.msg_in_new[pep as usize][node_neighbor_index] = self.msg_in[pep as usize][node_neighbor_index].clone();
            }
        }
    }
}