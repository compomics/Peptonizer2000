use std::collections::HashMap;
use minidom::Element;
use serde::Serialize;

#[derive(Debug, Serialize, Clone)]
pub struct Factor {
    // represents noisy OR cpds, has dimension n(parensports)xn(peptide states(=2))
    pub array: Vec<[f64; 2]>,
    pub array_labels: Vec<String>
}

#[derive(Debug, Serialize, Clone)]
pub enum Node {
    PeptideNode { id: String, initial_belief_0: f64, initial_belief_1: f64 },
    FactorNode { id: String, parent_number: i32, initial_belief: Factor },
    TaxonNode { id: String, initial_belief_0: f64, initial_belief_1: f64 },
    ConvolutionTreeNode { id: String, number_of_parents: i32 }
}

impl Node {
    pub fn id(&self) -> &str {
        match self {
            Node::PeptideNode { id, .. } => id,
            Node::FactorNode { id, .. } => id,
            Node::TaxonNode { id, .. } => id,
            Node::ConvolutionTreeNode { id, .. } => id
        }
    }

    pub fn category(&self) -> &str {
        match self {
            Node::PeptideNode { .. } => "peptide",
            Node::FactorNode { .. } => "factor",
            Node::TaxonNode { .. } => "taxon",
            Node::ConvolutionTreeNode { .. } => "convolution_tree"
        }
    }

    fn parse_data(data: &Element) -> (String, String) {
        let key = data.attr("key").unwrap().to_string();
        let value = data.text();
    
        (key, value)
    }

    pub fn parse_node(node: &Element) -> Result<Self, String> {
        // Process a node
        let id = node.attr("id").unwrap().to_string();
    
        // Initialize data for this node
        let mut current_node_data = HashMap::new();
        for data in node.children().filter(|d| d.name() == "data") {
            let (data_key, data_val) = Self::parse_data(data);
            current_node_data.insert(data_key, data_val);
        }
    
        match current_node_data.get("d2").map(String::as_str) {
            Some("factor") => {
                let parent_number: i32 = current_node_data.get("d3").unwrap().parse().unwrap();
                let initial_belief: Factor = Factor { array: Vec::new(), array_labels: Vec::new() };
                Ok(Node::FactorNode { id, parent_number, initial_belief })
            }
            Some("peptide") => {
                let initial_belief_0: f64 = current_node_data.get("d0").unwrap().parse().unwrap();
                let initial_belief_1: f64 = current_node_data.get("d1").unwrap().parse().unwrap();
                Ok(Node::PeptideNode { id, initial_belief_0, initial_belief_1 })
            }
            Some("taxon") => {
                let initial_belief_0: f64 = 0.0;
                let initial_belief_1: f64 = 0.0;
                Ok(Node::TaxonNode { id, initial_belief_0, initial_belief_1 })
            }
            _ => {
                Err("Node data has unknown type".to_string())
            }
        }
    }
}
