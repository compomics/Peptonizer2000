use petgraph::graph::{Graph, NodeIndex};
use petgraph::Undirected;
use crate::factor_graph_node::{Node};
use minidom::Element;
use std::collections::HashMap;

#[derive(Debug)]
pub struct CTFactorGraph {
    graph: Graph::<Node, (), Undirected>,
    node_map: HashMap<String, NodeIndex>
}

impl CTFactorGraph {
    pub fn to_string(&self) -> String {
        serde_json::to_string_pretty(&self.graph).unwrap()
    }

    fn parse_edge(edge: &Element) -> Result<(String, String), String> {
        let source: String = edge.attr("source").unwrap().to_string();
        let target: String = edge.attr("target").unwrap().to_string();
    
        Ok((source, target))
    }
    
    // Method to parse a GraphML string into the graph
    pub fn from_graphml(graphml_str: &str) -> Result<CTFactorGraph, String> {
        let mut graph = Graph::<Node, (), Undirected>::new_undirected();
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
    
                graph.add_edge(source_idx, target_idx, ());
            }
        }
    
        Ok( CTFactorGraph { graph, node_map })
    }
}
