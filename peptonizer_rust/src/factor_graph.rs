use petgraph::graph::{Graph, NodeIndex, EdgeIndex};
use petgraph::visit::{Dfs, EdgeRef};
use petgraph::Undirected;
use crate::factor_graph_node::{Factor, Node};
use minidom::Element;
use std::collections::{HashMap, HashSet};
use crate::utils::log;
use serde::Serialize;

#[derive(Debug, Serialize, Clone)]
pub struct EdgeData {
    message_length: Option<i32>
}

#[derive(Debug)]
pub struct CTFactorGraph {
    graph: Graph::<Node, EdgeData, Undirected>,
    node_map: HashMap<String, NodeIndex>
}

impl CTFactorGraph {
    pub fn to_string(&self) -> String {
        serde_json::to_string_pretty(&self.graph).unwrap()
    }

    pub fn add_node_categories(&self, node_categories: &mut HashMap<String, String>) {
        for node_index in self.graph.node_indices() {
            let node: &Node = self.graph.node_weight(node_index).unwrap();
            node_categories.insert(node.id().to_string(), node.category().to_string());
        }
    }

    pub fn node_count(&self) -> usize {
        self.graph.node_count()
    }

    fn parse_edge(edge: &Element) -> Result<(String, String), String> {
        let source: String = edge.attr("source").unwrap().to_string();
        let target: String = edge.attr("target").unwrap().to_string();
    
        Ok((source, target))
    }
    
    // Method to parse a GraphML string into the graph
    pub fn from_graphml(graphml_str: &str) -> Result<CTFactorGraph, String> {
        let mut graph = Graph::<Node, EdgeData, Undirected>::new_undirected();
        let mut node_map: HashMap<String, NodeIndex> = HashMap::new();
    
        let root: Element = graphml_str.parse().unwrap();
    
        for graph_xml in root.children().filter(|n| n.name() == "graph") {
            for node_xml in graph_xml.children().filter(|n| n.name() == "node") {
                let node: Node = Node::parse_node(node_xml).unwrap();
    
                let node_id: String = node.id().to_string();
                let node_index: NodeIndex = graph.add_node(node);
                node_map.insert(node_id, node_index);
            }
    
            for edge_xml in graph_xml.children().filter(|n| n.name() == "edge") {
                let (source, target) = Self::parse_edge(edge_xml).unwrap();
    
                let source_idx: NodeIndex = *node_map.get(&source).unwrap();
                let target_idx: NodeIndex = *node_map.get(&target).unwrap();
                let edge_data: EdgeData = EdgeData { message_length: None };
    
                graph.add_edge(source_idx, target_idx, edge_data);
            }
        }
    
        Ok( CTFactorGraph { graph, node_map })
    }

    pub fn fill_in_priors(&mut self, prior: f64) {
        for node_index in self.graph.node_indices() {
            let node: &mut Node = self.graph.node_weight_mut(node_index).unwrap();

            if let Node::TaxonNode { id, .. } = &node {
                *node = Node::TaxonNode { id: id.to_string(), initial_belief_0: 1.0 - prior, initial_belief_1: prior };
            }
        }
    }

    pub fn fill_in_factors(&mut self, alpha: f64, beta: f64, regularized: bool) {
        for node_index in self.graph.node_indices() {
            let node: &mut Node = self.graph.node_weight_mut(node_index).unwrap();

            if let Node::FactorNode { id, parent_number, .. } = &node {
                let degree: i32 = *parent_number;

                let mut cpd_array: Vec<[f64; 2]> = Vec::with_capacity(degree as usize + 1);
                let mut cpd_array_regularized = cpd_array.clone();
                let exponent_array: Vec<i32> = (0..=degree).collect();
                let divide_array: Vec<f64> = std::iter::once(1i32).chain(1..=degree).map(|x| x as f64).collect();
                
                // regularize cpd priors to penalize higher number of parents
                // log domain to avoid underflow
                let mut cpd_sum: f64 = 0.0;
                let mut cpd_regularized_sum: f64 = 0.0;
                for (i, exp) in exponent_array.iter().enumerate() {
                    let cpd_0 = (1.0 - alpha).powi(*exp) * (1.0 - beta);
                    let cpd_1 = 1.0 - cpd_0;
                    cpd_sum += cpd_0 + cpd_1;
                    cpd_array.push([cpd_0, cpd_1]);

                    let cpd_regularized_0 = (cpd_0.powi(*exp) * (1.0 - beta)) / divide_array[i];
                    let cpd_regularized_1 = 1.0 - cpd_regularized_0;
                    cpd_regularized_sum += cpd_regularized_0 + cpd_regularized_1;
                    cpd_array_regularized.push([cpd_regularized_0, cpd_regularized_1]);
                }

                // Normalize arrays (assuming normalize and avoid_underflow are implemented)
                Self::normalize_cpd(&mut cpd_array, cpd_sum, false);
                Self::normalize_cpd(&mut cpd_array_regularized, cpd_regularized_sum, true);
                
                // Create factor
                let factor_to_add = if regularized {
                    Factor {
                        array: cpd_array_regularized,
                        array_labels: vec![format!("placeholder"), format!("{}0", node_index.index()), format!("{}1", node_index.index())],
                    }
                } else {
                    Factor {
                        array: cpd_array,
                        array_labels : vec![format!("placeholder"), format!("{}0", node_index.index()), format!("{}1", node_index.index())],
                    }
                };
                
                // Add factor to the node's attributes
                *node = Node::FactorNode { id: id.to_string(), parent_number: *parent_number, initial_belief: factor_to_add };
            }
        }
    }

    fn normalize_cpd(arr: &mut Vec<[f64; 2]>, sum: f64, avoid_underflow: bool) {
        for cpd in arr.iter_mut() {
            cpd[0] /= sum;
            cpd[1] /= sum;
    
            if avoid_underflow {
                if cpd[0] < 1e-30 {
                    cpd[0] = 1e-30;
                }
                if cpd[1] < 1e-30 {
                    cpd[1] = 1e-30;
                }
            }
        }
    }

    pub fn add_ct_nodes(&mut self) {
        // When creating the CTGraph and not just reading from a previously saved graph format, use this function to add the CT nodes
        
        let mut edges_to_add: Vec<(NodeIndex, NodeIndex, EdgeData)> = Vec::new();
        let mut edges_to_remove: Vec<(NodeIndex, NodeIndex)> = Vec::new();

        for node_index in self.graph.node_indices() {
            let node: &Node = self.graph.node_weight(node_index).unwrap();

            if let Node::FactorNode { .. } = &node {
                if self.graph.neighbors(node_index).count() > 2 {

                    let mut prot_ids: Vec<String> = Vec::new();
                    let mut prot_list: Vec<NodeIndex> = Vec::new();
                    for neighbor_index in self.graph.neighbors(node_index) {

                        let neighbor: &Node = self.graph.node_weight(neighbor_index).unwrap();
                        if let Node::FactorNode { id, .. } = &neighbor {
                            prot_list.push(neighbor_index);
                            prot_ids.push(id.to_string());
                        }
                    }

                    let node_id = prot_ids.join(" ");
                    let new_node = Node::ConvolutionTreeNode { id: node_id.clone(), number_of_parents: prot_list.len() as i32 };
                    let new_node_index: NodeIndex = self.graph.add_node(new_node);
                    self.node_map.insert(node_id, new_node_index);

                    let edge_data = EdgeData { message_length: Some(prot_list.len() as i32 + 1) };
                    edges_to_add.push((new_node_index, node_index, edge_data));
                    for neighbor in prot_list {
                        let edge_data = EdgeData { message_length: None };
                        edges_to_add.push((new_node_index, neighbor, edge_data));
                        edges_to_remove.push((node_index, neighbor));
                    }
                }
                
            }
        }

        for (node1, node2, edge_data) in edges_to_add {
            self.graph.add_edge(node1, node2, edge_data);
        }

        for (node1, node2) in edges_to_remove  {
            let edge_id: EdgeIndex = self.graph.find_edge(node1, node2).unwrap();
            self.graph.remove_edge(edge_id);
        }
    }

    /// Finds the connected components in an undirected graph and returns a Vec of Vecs containing nodes in each component.
    pub fn connected_components(&self) -> Vec<Self> {
        let mut visited: HashSet<NodeIndex> = HashSet::new();
        let mut components: Vec<Self> = Vec::new();

        for start_node in self.graph.node_indices() {
            if !visited.contains(&start_node) {
                let mut component_nodes: HashSet<NodeIndex> = HashSet::new();
                let mut old_to_new_nodes: HashMap<NodeIndex, NodeIndex> = HashMap::new();

                let mut subgraph = Graph::<Node, EdgeData, Undirected>::new_undirected();
                let mut node_map: HashMap<String, NodeIndex> = HashMap::new();

                let node_weight: Node = self.graph[start_node].clone();
                let new_node: NodeIndex = subgraph.add_node(node_weight.clone());
                node_map.insert(node_weight.id().to_string(), new_node);

                let mut dfs = Dfs::new(&self.graph, start_node);
                while let Some(node) = dfs.next(&self.graph) {
                    if visited.insert(node) {
                        component_nodes.insert(node);

                        let node_weight: Node = self.graph[node].clone();
                        let new_node: NodeIndex = subgraph.add_node(node_weight.clone());
                        node_map.insert(node_weight.id().to_string(), new_node);

                        old_to_new_nodes.insert(node, new_node);
                    }
                }

                for edge in self.graph.edge_references() {
                    let (source, target): (NodeIndex, NodeIndex) = (edge.source(), edge.target());
                    if component_nodes.contains(&source) && component_nodes.contains(&target) {
                        let weight = edge.weight().clone();
                        subgraph.add_edge(old_to_new_nodes[&source], old_to_new_nodes[&target], weight);
                    }
                }

                let factor_subgraph = Self { graph: subgraph, node_map };
                components.push(factor_subgraph);
            }
        }

        components
    }
}
