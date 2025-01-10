use crate::node::{Factor, Node};
use minidom::Element;
use std::collections::{HashMap, HashSet};
use crate::utils::log;
use serde::Serialize;

#[derive(Debug, Serialize, Clone)]
pub struct Edge {
    id: i32,
    node1_id: i32,
    node2_id: i32,
    message_length: Option<i32>
}

impl Edge {

    pub fn get_id(&self) -> i32 {
        self.id
    }

    pub fn get_nodes(&self) -> (i32, i32) {
        (self.node1_id, self.node2_id)
    }

    pub fn get_message_length(&self) -> Option<i32> {
        self.message_length
    }
}

#[derive(Debug)]
pub struct CTFactorGraph {
    nodes: Vec<Node>,
    edges: Vec<Edge>,
}

impl CTFactorGraph {

    pub fn add_node_categories(&self, node_categories: &mut HashMap<String, String>) {
        for node in &self.nodes {
            node_categories.insert(node.get_name().to_string(), node.category().to_string());
        }
    } 

    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    pub fn edge_count(&self) -> usize {
        self.edges.len()
    }

    fn parse_edge(edge: &Element) -> Result<(String, String), String> {
        let source: String = edge.attr("source").unwrap().to_string();
        let target: String = edge.attr("target").unwrap().to_string();
    
        Ok((source, target))
    }
    
    // Method to parse a GraphML string into the graph
    pub fn from_graphml(graphml_str: &str) -> Result<CTFactorGraph, String> {
        let mut nodes: Vec<Node> = Vec::new();
        let mut edges: Vec<Edge> = Vec::new();
        let mut node_map: HashMap<String, i32> = HashMap::new();
    
        let root: Element = graphml_str.parse().unwrap();
        
        let mut next_node_id = 0;
        let mut next_edge_id = 0;
        for graph_xml in root.children().filter(|n| n.name() == "graph") {
            for node_xml in graph_xml.children().filter(|n| n.name() == "node") {
                let node: Node = Node::parse_node(node_xml, next_node_id).unwrap();
                let node_name: String = node.get_name().to_string();
                node_map.insert(node_name, next_node_id);
                next_node_id += 1;

                nodes.push(node);
            }
    
            for edge_xml in graph_xml.children().filter(|n| n.name() == "edge") {
                let (source, target) = Self::parse_edge(edge_xml).unwrap();
    
                let node1_id: i32 = *node_map.get(&source).unwrap();
                let node2_id: i32 = *node_map.get(&target).unwrap();
                let edge: Edge = Edge { id: next_edge_id, node1_id, node2_id, message_length: None };
                next_edge_id += 1;
    
                let node1: &mut Node = &mut nodes[node1_id as usize];
                node1.add_incident_edge(edge.get_id());
                let node2: &mut Node = &mut nodes[node2_id as usize];
                node2.add_incident_edge(edge.get_id());
                edges.push(edge);
            }
        }
    
        Ok( CTFactorGraph { nodes, edges })
    }

    pub fn fill_in_priors(&mut self, prior: f64) {
        for node in &mut self.nodes {
            node.fill_in_prior(prior);
        }
    }

    pub fn fill_in_factors(&mut self, alpha: f64, beta: f64, regularized: bool) {
        for node in &mut self.nodes {
            node.fill_in_factor(alpha, beta, regularized);
        }
    }

    pub fn get_neighbors(&self, node: &Node) -> Vec<i32> {
        let mut neighbors = Vec::with_capacity(node.neighbors_count() as usize);
        for edge_id in node.get_incident_edges() {
            let (node1_id, node2_id) = self.edges[*edge_id as usize].get_nodes();
            let neighbor: i32 = if node1_id == node.get_id() { node2_id } else { node1_id };
            neighbors.push(neighbor);
        }
        
        neighbors
    }

    pub fn add_ct_nodes(&mut self) {
        // When creating the CTGraph and not just reading from a previously saved graph format, use this function to add the CT nodes
        
        let mut edges_to_add: Vec<Edge> = Vec::new();
        let mut edges_to_remove: HashSet<(i32, i32)> = HashSet::new();
        let mut nodes_to_add: Vec<Node> = Vec::new();

        // Add nodes and keep track of edges to add/remove
        let mut next_edge_id: i32 = self.edges.len() as i32;
        let mut next_node_id: i32 = self.nodes.len() as i32;
        for node in &self.nodes {
            if node.is_factor_node() {
                if node.neighbors_count() > 2 {

                    let mut prot_names: Vec<String> = Vec::new();
                    let mut prot_ids: Vec<i32> = Vec::new();
                    for neighbor_id in self.get_neighbors(node) {
                        let neighbor: &Node = &self.nodes[neighbor_id as usize];
                        if neighbor.is_factor_node() {
                            prot_ids.push(neighbor_id);
                            prot_names.push(neighbor.get_name().to_string());
                        }
                    }

                    let new_node_name = prot_names.join(" ");
                    let new_node_id = next_node_id;
                    let new_node = Node::new_convolution_node(new_node_id, new_node_name, prot_ids.len() as i32);
                    next_node_id += 1;
                    nodes_to_add.push(new_node);

                    let edge = Edge { id: next_edge_id, node1_id: new_node_id, node2_id: node.get_id(), message_length: Some(prot_ids.len() as i32 + 1) };
                    next_edge_id += 1;
                    edges_to_add.push(edge);
                    for neighbor_id in prot_ids {
                        let edge = Edge { id: next_edge_id, node1_id: new_node_id, node2_id: neighbor_id, message_length: None };
                        next_edge_id += 1;
                        edges_to_add.push(edge);
                        edges_to_remove.insert((node.get_id(), neighbor_id));
                        edges_to_remove.insert((neighbor_id, node.get_id()));
                    }
                }
                
            }
        }

        // Remove edges
        let mut new_edges: Vec<Edge> = Vec::with_capacity(self.edges.len() + edges_to_add.len() - edges_to_remove.len());
        let mut next_edge_id = 0;
        for edge in &self.edges {
            if ! edges_to_remove.contains(&(edge.node1_id, edge.node2_id)) {
                let mut new_edge = edge.clone();
                new_edge.id = next_edge_id;
                next_edge_id += 1;
                new_edges.push(new_edge);
            }
        }
        // Add new edges
        for edge in edges_to_add {
            let mut new_edge = edge.clone();
            new_edge.id = next_edge_id;
            next_edge_id += 1;
            new_edges.push(new_edge);
        }

        // Update the incident edges in the nodes
        let mut new_nodes: Vec<Node> = Vec::with_capacity(self.nodes.len() + nodes_to_add.len());
        for node in &self.nodes {
            new_nodes.push(Node::new(node.get_id(), node.get_name().to_string(), node.get_subtype().clone()));
        }
        for node in nodes_to_add {
            new_nodes.push(node);
        }
        
        for edge in &new_edges {
            let (node1_id, node2_id) = edge.get_nodes();
            new_nodes[node1_id as usize].add_incident_edge(edge.get_id());
            new_nodes[node2_id as usize].add_incident_edge(edge.get_id());
        }

        self.nodes = new_nodes;
        self.edges = new_edges;
    }

    /// Finds the connected components in an undirected graph and returns a Vec of Vecs containing nodes in each component.
    pub fn connected_components(&self) -> Vec<Self> {
        let mut visited: HashSet<i32> = HashSet::new();
        let mut components: Vec<Self> = Vec::new();

        for start_node in &self.nodes {
            if visited.insert(start_node.get_id()) {
                let mut component_ids: Vec<i32> = Vec::new();
                let mut old_to_new_nodes: HashMap<i32, i32> = HashMap::new();

                let mut new_nodes: Vec<Node> = Vec::new();
                let mut new_edges: Vec<Edge> = Vec::new();

                // Find ids of nodes to include in component
                component_ids.push(start_node.get_id());
                old_to_new_nodes.insert(start_node.get_id(), 0);
                self.find_component_rec(start_node.get_id(), &mut component_ids, &mut old_to_new_nodes, &mut visited);

                // log(&format!("{}", component_ids.len())); = 24000

                // Create new nodes
                for node_id in &component_ids {
                    let node = self.nodes[*node_id as usize].copy_with_id(old_to_new_nodes[&node_id]);
                    new_nodes.push(node);
                }

                // Select edges to keep and update the node ids
                let mut next_edge_id: i32 = 0;
                let mut component_edge_ids: HashSet<i32> = HashSet::new();
                let mut old_to_new_edges: HashMap<i32, i32> = HashMap::new();
                for edge in &self.edges {

                    let (source, target): (i32, i32) = edge.get_nodes();
                    if component_ids.contains(&source) && component_ids.contains(&target) {

                        let (new_source, new_target): (i32, i32) = (old_to_new_nodes[&source], old_to_new_nodes[&target]);
                        let new_edge = Edge { id: next_edge_id, node1_id: new_source, node2_id: new_target, message_length: edge.get_message_length() };
                        next_edge_id += 1;

                        component_edge_ids.insert(edge.get_id());
                        old_to_new_edges.insert(edge.get_id(), new_edge.get_id());

                        new_edges.push(new_edge);
                    }
                }

                // Update edge ids of incident edges
                for node in &mut new_nodes {
                    let new_incident_edges: Vec<i32> = node.get_incident_edges().into_iter().filter(|e| component_edge_ids.contains(e)).map(|e| old_to_new_edges[e]).collect();
                    node.set_incident_edges(new_incident_edges);
                }
                
                log(&format!("{}", new_nodes.len()));
                // Create graph and add to components
                let subgraph = Self { nodes: new_nodes, edges: new_edges };
                components.push(subgraph);
            }
        }

        components
    }

    fn find_component_rec(
        &self, 
        start_id: i32, 
        component_ids: &mut Vec<i32>, 
        old_to_new_nodes: &mut HashMap<i32, i32>, 
        visited: &mut HashSet<i32>
    ) {
        let start_node: &Node = &self.nodes[start_id as usize];
        for neighbor_id in self.get_neighbors(&start_node) {
            if visited.insert(neighbor_id) {
                let next_id: i32 = component_ids.len() as i32;
                component_ids.push(neighbor_id);
                old_to_new_nodes.insert(neighbor_id, next_id);
                self.find_component_rec(neighbor_id, component_ids, old_to_new_nodes, visited);                
            }
        }
    }
}
