//! The main FluidNetwork implementation.
//!
//! Uses Vec-based storage with custom indices for safety.

use crate::elements::{CurrentNodeState, Edge, EdgeType, ExternalInputs, FlowResult, Node};
use crate::indices::{EdgeIndex, NodeIndex};
use std::collections::HashSet;
use std::f32::consts::E;

/// The core fluid network structure.
#[derive(Debug)]
pub struct FluidNetwork {
    pub nodes: Vec<Node>,
    pub edges: Vec<Edge>,
    /// For each node, list of connected edge indices.
    pub node_edges: Vec<Vec<EdgeIndex>>,
    /// For each edge, the connected node indices.
    pub edge_nodes: Vec<(NodeIndex, NodeIndex)>,
    /// Current runtime state of each node.
    pub current_node_states: Vec<CurrentNodeState>,
    /// Current runtime flow of each edge.
    pub current_flows: Vec<FlowResult>,
    /// Used names for uniqueness.
    pub used_names: HashSet<String>,
    /// Counter for auto-generated names.
    name_counter: u32,
}

impl FluidNetwork {
    /// Creates a new empty FluidNetwork.
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            edges: Vec::new(),
            node_edges: Vec::new(),
            edge_nodes: Vec::new(),
            current_node_states: Vec::new(),
            current_flows: Vec::new(),
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

    /// Initializes the runtime state for a node.
    fn initialize_node_state(&self, node: &Node) -> CurrentNodeState {
        match node {
            Node::Thermal(node) => CurrentNodeState::Thermal {
                temperature: node.initial_temp,
            },
            Node::Fluid(node) => CurrentNodeState::Fluid {
                liquid_temperature: None, // Initial as None, can be set later
                gas_temperature: None,
                pressure: 101325.0, // Atmospheric, placeholder
                volume: node.volume,
                liquid_fraction: node.liquid_fraction,
            },
            Node::ThermalBoundary(node) => CurrentNodeState::ThermalBoundary {
                temperature: node.temperature,
            },
            Node::FluidBoundary(node) => CurrentNodeState::FluidBoundary {
                liquid_temperature: Some(node.temperature),
                gas_temperature: None, // Placeholder
                pressure: node.pressure,
                liquid_fraction: 0.,
            },
            Node::FluidJunction(_node) => CurrentNodeState::Fluid {
                liquid_temperature: None,
                gas_temperature: None,
                pressure: 101325.0,   // Placeholder
                volume: 0.0,          // To be updated
                liquid_fraction: 0.5, // Placeholder
            },
        }
    }

    /// Adds a node and returns its index.
    pub fn add_node(&mut self, mut node: Node, name: Option<&str>) -> NodeIndex {
        // Set the name on the node
        self.set_node_name(&mut node, name);
        // Initialize state
        let initial_state = self.initialize_node_state(&node);
        self.nodes.push(node);
        self.current_node_states.push(initial_state);
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
        let edge_type = edge.edge_type.clone();
        // Set name on edge
        edge.name = Some(self.generate_unique_name(name, "Edge"));
        self.edges.push(edge);
        self.edge_nodes.push((node1, node2));
        let edge_index = EdgeIndex::new(self.edges.len() - 1, self.edges.len()).unwrap();

        // Add to node connections
        self.node_edges[node1.index()].push(edge_index);
        self.node_edges[node2.index()].push(edge_index);

        match edge_type {
            EdgeType::Fluid => self.current_flows.push(FlowResult::Fluid {
                liquid_mass_flow: 0.,
                gas_mass_flow: 0.,
                liquid_specific_enthalpy: 4200., // Dummy, should be based on fluid props
                gas_specific_enthalpy: 2000.,    // Dummy, should be based on fluid props
            }),
            EdgeType::Thermal => self
                .current_flows
                .push(FlowResult::Thermal { thermal_power: 0. }),
        }
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

    /// Phases the network by one time step with given external inputs and time step.
    pub fn step(&mut self, external_inputs: &ExternalInputs, dt: f64) {
        // Collect all flow results into `self.current_flows`
        for (edge_idx, edge) in self.edges.iter().enumerate() {
            let edge_idx_checked = EdgeIndex::new(edge_idx, self.edges.len()).unwrap();
            let (node_a_idx, node_b_idx) = self.edge_nodes[edge_idx];
            let node_a_state = &self.current_node_states[node_a_idx.index()];
            let node_b_state = &self.current_node_states[node_b_idx.index()];
            let flow_result = edge.behavior.compute_flow(
                node_a_idx,
                node_a_state,
                node_b_idx,
                node_b_state,
                external_inputs,
            );
            self.current_flows[edge_idx_checked.index()] = flow_result;
        }

        // Update node states (simple integration for thermal nodes)
        for (node_idx, node_state) in self.current_node_states.iter_mut().enumerate() {
            let node_idx_checked = NodeIndex::new(node_idx, self.nodes.len()).unwrap();
            if let CurrentNodeState::Thermal { temperature } = node_state {
                if let Some(node) = self.nodes.get(node_idx) {
                    if let Node::Thermal(thermal_node) = node {
                        let mut net_power = 0.0;
                        // Sum incoming/outgoing thermal power from connected edges
                        for &edge_idx in &self.node_edges[node_idx] {
                            let edge_idx = edge_idx.index();
                            let flow = &self.current_flows[edge_idx];
                            if let FlowResult::Thermal { thermal_power } = flow {
                                let (n1, n2) = self.edge_nodes[edge_idx];
                                // If this node is A, power is positive A to B, else negative
                                if n1 == node_idx_checked {
                                    net_power += thermal_power;
                                } else {
                                    net_power -= thermal_power;
                                }
                            }
                        }
                        // Update temperature: dT = power * dt / mass
                        *temperature += net_power * dt / thermal_node.thermal_mass;
                    }
                }
            }
        }
    }

    /// Runs the simulation for a number of steps with constant external inputs and dt.
    pub fn run(&mut self, external_inputs: &ExternalInputs, dt: f64, num_steps: usize) {
        for _ in 0..num_steps {
            self.step(external_inputs, dt);
        }
    }

    // TODO: Implement simplification and validation methods.
}
