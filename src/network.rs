//! The main FluidNetwork implementation.
//!
//! Uses Vec-based storage with custom indices for safety.

use crate::elements::{Edge, Node};
use crate::indices::{EdgeIndex, NodeIndex};
use std::collections::HashSet;

/// The core fluid network structure.
#[derive(Debug)]
pub struct FluidNetwork {
    pub nodes: Vec<Node>,
    pub edges: Vec<Edge>,
    /// For each node, list of connected edge indices.
    pub node_edges: Vec<Vec<EdgeIndex>>,
    /// For each edge, the connected node indices.
    pub edge_nodes: Vec<(NodeIndex, NodeIndex)>,
    /// Used names for uniqueness.
    pub used_names: HashSet<String>,
    /// Counter for auto-generated names.
    pub name_counter: u32,
}

impl FluidNetwork {
    /// Creates a new empty FluidNetwork.
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            edges: Vec::new(),
            node_edges: Vec::new(),
            edge_nodes: Vec::new(),
            used_names: HashSet::new(),
            name_counter: 1,
        }
    }

    /// Generates a unique name for an element.
    /// If proposed_name is provided, uses it if unique, else adds suffix.
    /// If None, generates a default name based on element_type.
    fn generate_unique_name(&mut self, proposed_name: Option<&str>, element_type: &str) -> String {
        if let Some(prop) = proposed_name {
            if !self.used_names.contains(prop) {
                self.used_names.insert(prop.to_string());
                return prop.to_string();
            } else {
                let mut suffix = 1;
                loop {
                    let new_name = format!("{}_{}", prop, suffix);
                    if !self.used_names.contains(&new_name) {
                        self.used_names.insert(new_name.clone());
                        return new_name;
                    }
                    suffix += 1;
                }
            }
        } else {
            loop {
                let name = format!("{}_{}", element_type, self.name_counter);
                self.name_counter += 1;
                if !self.used_names.contains(&name) {
                    self.used_names.insert(name.clone());
                    return name;
                }
            }
        }
    }

    /// Sets the name on a Node.
    fn set_node_name(&mut self, node: &mut Node, name: Option<&str>) {
        match node {
            Node::Thermal(n) => {
                n.name = Some(self.generate_unique_name(name, "Thermal"));
            }
            Node::Fluid(n) => {
                n.name = Some(self.generate_unique_name(name, "Fluid"));
            }
            Node::ThermalBoundary(n) => {
                n.name = Some(self.generate_unique_name(name, "ThermalBoundary"));
            }
            Node::FluidBoundary(n) => {
                n.name = Some(self.generate_unique_name(name, "FluidBoundary"));
            }
            Node::FluidJunction(n) => {
                n.name = Some(self.generate_unique_name(name, "FluidJunction"));
            }
        }
    }

    /// Adds a node and returns its index.
    pub fn add_node(&mut self, mut node: Node, name: Option<&str>) -> NodeIndex {
        // Set the name on the node
        self.set_node_name(&mut node, name);
        self.nodes.push(node);
        self.node_edges.push(Vec::new()); // Empty connections initially
        NodeIndex::new(self.nodes.len() - 1, self.nodes.len()).unwrap()
    }

    /// Adds an edge between two nodes and returns its index.
    pub fn add_edge(
        &mut self,
        mut edge: Edge,
        node1: NodeIndex,
        node2: NodeIndex,
        name: Option<&str>,
    ) -> Option<EdgeIndex> {
        if node1.index() >= self.nodes.len() || node2.index() >= self.nodes.len() {
            return None; // Invalid indices
        }
        // Set name on edge
        edge.name = Some(self.generate_unique_name(name, "Edge"));
        self.edges.push(edge);
        self.edge_nodes.push((node1, node2));
        let edge_index = EdgeIndex::new(self.edges.len() - 1, self.edges.len()).unwrap();

        // Add to node connections
        self.node_edges[node1.index()].push(edge_index);
        self.node_edges[node2.index()].push(edge_index);

        Some(edge_index)
    }

    /// Looks up a node by name and returns its index if found.
    pub fn find_node_by_name(&self, name: &str) -> Option<NodeIndex> {
        for (i, node) in self.nodes.iter().enumerate() {
            let node_name = match node {
                Node::Thermal(n) => &n.name,
                Node::Fluid(n) => &n.name,
                Node::ThermalBoundary(n) => &n.name,
                Node::FluidBoundary(n) => &n.name,
                Node::FluidJunction(n) => &n.name,
            };
            if let Some(n) = node_name {
                if n == name {
                    return NodeIndex::new(i, self.nodes.len());
                }
            }
        }
        None
    }

    /// Looks up an edge by name and returns its index if found.
    pub fn find_edge_by_name(&self, name: &str) -> Option<EdgeIndex> {
        for (i, edge) in self.edges.iter().enumerate() {
            if let Some(edge_name) = &edge.name {
                if edge_name == name {
                    return EdgeIndex::new(i, self.edges.len());
                }
            }
        }
        None
    }

    // TODO: Implement simplification and validation methods.
    // TODO: Implement execution methods: step, run, etc.
    // These will be added in subsequent phases.
}
