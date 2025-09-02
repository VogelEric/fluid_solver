//! Defines the core elements of the fluid network.
//!
//! This module contains structs and enums for nodes, edges, and their properties.

use crate::indices::{EdgeIndex, NodeIndex};
use std::collections::HashMap;

mod hydrostatic_tank;
pub use hydrostatic_tank::HydrostaticTank;

#[derive(Debug, Clone)]
pub struct ThermalNode {
    pub name: Option<String>,
    pub thermal_mass: f64,
    pub initial_temp: f64,
}

/// Represents a fluid node with internal state.
#[derive(Debug, Clone)]
pub struct FluidNode {
    pub name: Option<String>,
    pub volume: f64,
    pub liquid_fraction: f64,
}

/// Represents a thermal boundary condition.
#[derive(Debug, Clone)]
pub struct ThermalBoundary {
    pub name: Option<String>,
    pub temperature: f64, // External input
}

/// Represents a fluid boundary condition.
#[derive(Debug, Clone)]
pub struct FluidBoundary {
    pub name: Option<String>,
    pub temperature: f64,
    pub pressure: f64,
}

/// Represents a placeholder for complex networks.
#[derive(Debug, Clone)]
pub struct FluidJunction {
    pub name: Option<String>,
    pub estimated_volume: f64, // Placeholder, replaced during simplification
}

/// Enum for all node types.
#[derive(Debug, Clone)]
pub enum Node {
    Thermal(ThermalNode),
    Fluid(FluidNode),
    ThermalBoundary(ThermalBoundary),
    FluidBoundary(FluidBoundary),
    FluidJunction(FluidJunction),
    HydrostaticTank(HydrostaticTank),
}

impl Node {
    pub fn compute_new_node_state(
        &self,
        current_state: &CurrentNodeState,
        connected_flows: &Vec<(EdgeIndex, FlowResult)>,
        external_inputs: &ExternalInputs,
        node_index: NodeIndex,
        dt: f64,
    ) -> CurrentNodeState {
        match self {
            Node::Thermal(node) => node.compute_new_state(
                current_state,
                connected_flows,
                external_inputs,
                node_index,
                dt,
            ),
            Node::Fluid(node) => node.compute_new_state(
                current_state,
                connected_flows,
                external_inputs,
                node_index,
                dt,
            ),
            Node::ThermalBoundary(node) => node.compute_new_state(
                current_state,
                connected_flows,
                external_inputs,
                node_index,
                dt,
            ),
            Node::FluidBoundary(node) => node.compute_new_state(
                current_state,
                connected_flows,
                external_inputs,
                node_index,
                dt,
            ),
            Node::FluidJunction(node) => node.compute_new_state(
                current_state,
                connected_flows,
                external_inputs,
                node_index,
                dt,
            ),
            Node::HydrostaticTank(node) => node.compute_new_state(
                current_state,
                connected_flows,
                external_inputs,
                node_index,
                dt,
            ),
        }
    }
}

/// Represents external inputs that can override node states at runtime.
#[derive(Debug, Clone)]
pub struct ExternalInputs {
    pub temperature_overrides: HashMap<NodeIndex, f64>,
    pub pressure_overrides: HashMap<NodeIndex, f64>,
}

/// Represents the current runtime state of a node.
/// Separate temperatures for liquid and gas phases in fluid nodes.
#[derive(Debug, Clone, PartialEq)]
pub enum CurrentNodeState {
    Thermal {
        temperature: f64,
    },
    Fluid {
        liquid_temperature: Option<f64>,
        gas_temperature: Option<f64>,
        pressure: f64,
        volume: f64,
        liquid_fraction: f64,
    },
    ThermalBoundary {
        temperature: f64,
    },
    FluidBoundary {
        liquid_temperature: Option<f64>,
        gas_temperature: Option<f64>,
        pressure: f64,
        liquid_fraction: f64,
    },
    HydrostaticTank {
        liquid_volume: f64,
        liquid_height: f64,
        gas_pressure_surface: f64,
        temperature: f64,
    },
}

/// Result of edge flow computation.
#[derive(Debug, Clone)]
pub enum FlowResult {
    Thermal {
        thermal_power: f64,
    },
    Fluid {
        liquid_mass_flow: f64,
        gas_mass_flow: f64,
        liquid_specific_enthalpy: f64,
        gas_specific_enthalpy: f64,
    },
}

/// High level edge classification
#[derive(Debug, Clone)]
pub enum EdgeType {
    Thermal,
    Fluid,
}

/// Refactored Edge struct that uses trait objects for flexibility.
pub struct Edge {
    pub name: Option<String>,
    pub behavior: Box<dyn EdgeBehavior>,
    pub edge_type: EdgeType,
    /// Heights where this edge connects to its nodes (relative to node reference point)
    pub connection_height_a: f64,
    pub connection_height_b: f64,
}

impl std::fmt::Debug for Edge {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Edge {{ ... }}")
    }
}

/// Trait for edge behavior
pub trait EdgeBehavior {
    fn compute_flow(
        &self,
        node_a_index: NodeIndex,
        node_a_state: &CurrentNodeState,
        connection_height_a: f64,
        node_b_index: NodeIndex,
        node_b_state: &CurrentNodeState,
        connection_height_b: f64,
        external_inputs: &ExternalInputs,
    ) -> FlowResult;

    /// Returns the temperature at nodes a and b at either end of this node
    /// format is (temp_a, temp_b)
    fn get_connected_edge_temps(
        &self,
        node_a_state: &CurrentNodeState,
        node_b_state: &CurrentNodeState,
    ) -> (f64, f64) {
        let temp_a = match node_a_state {
            CurrentNodeState::Thermal { temperature } => *temperature,
            CurrentNodeState::ThermalBoundary { temperature } => *temperature,
            CurrentNodeState::Fluid {
                liquid_temperature,
                gas_temperature,
                liquid_fraction,
                ..
            } => match (*liquid_temperature, *gas_temperature) {
                (Some(l_temp), None) => l_temp,
                (None, Some(g_temp)) => g_temp,
                (Some(l_temp), Some(g_temp)) => {
                    g_temp * (1. - liquid_fraction) + l_temp * (liquid_fraction)
                }
                (None, None) => unreachable!(
                    "A Fluid state should always have a temperature for one of the phases"
                ),
            },
            _ => {
                unreachable!("Validation should ensure all edge connecitons are valid.");
            } // Invalid, should not happen for thermal edges
        };
        let temp_b = match node_b_state {
            CurrentNodeState::Thermal { temperature } => *temperature,
            CurrentNodeState::ThermalBoundary { temperature } => *temperature,
            CurrentNodeState::Fluid {
                liquid_temperature,
                gas_temperature,
                liquid_fraction,
                ..
            } => match (*liquid_temperature, *gas_temperature) {
                (Some(l_temp), None) => l_temp,
                (None, Some(g_temp)) => g_temp,
                (Some(l_temp), Some(g_temp)) => {
                    g_temp * (1. - liquid_fraction) + l_temp * (liquid_fraction)
                }
                (None, None) => unreachable!(
                    "A Fluid state should always have a temperature for one of the phases"
                ),
            },
            _ => {
                unreachable!("Validation should ensure all edge connecitons are valid.");
            } // Invalid, should not happen for thermal edges
        };
        (temp_a, temp_b)
    }

    fn get_connected_fluid_states(
        &self,
        node_a_state: &CurrentNodeState,
        node_b_state: &CurrentNodeState,
    ) -> (
        (f64, Option<f64>, Option<f64>, f64),
        (f64, Option<f64>, Option<f64>, f64),
    ) {
        let (pressure_a, temp_a_liquid, temp_a_gas, liquid_fraction_a) = match node_a_state {
            CurrentNodeState::Fluid {
                liquid_temperature,
                gas_temperature,
                pressure,
                liquid_fraction,
                ..
            } => (
                *pressure,
                *liquid_temperature,
                *gas_temperature,
                *liquid_fraction,
            ),
            CurrentNodeState::FluidBoundary {
                liquid_temperature,
                gas_temperature,
                pressure,
                liquid_fraction,
            } => (
                *pressure,
                *liquid_temperature,
                *gas_temperature,
                *liquid_fraction,
            ),
            _ => {
                unreachable!("Validation should ensure all edge connecitons are valid.");
            } // Invalid for fluid edges
        };
        let (pressure_b, temp_b_liquid, temp_b_gas, liquid_fraction_b) = match node_b_state {
            CurrentNodeState::Fluid {
                liquid_temperature,
                gas_temperature,
                pressure,
                liquid_fraction,
                ..
            } => (
                *pressure,
                *liquid_temperature,
                *gas_temperature,
                *liquid_fraction,
            ),
            CurrentNodeState::FluidBoundary {
                liquid_temperature,
                gas_temperature,
                pressure,
                liquid_fraction,
            } => (
                *pressure,
                *liquid_temperature,
                *gas_temperature,
                *liquid_fraction,
            ),
            _ => {
                unreachable!("Validation should ensure all edge connecitons are valid.");
            } // Invalid for fluid edges
        };
        (
            (pressure_a, temp_a_liquid, temp_a_gas, liquid_fraction_a),
            (pressure_b, temp_b_liquid, temp_b_gas, liquid_fraction_b),
        )
    }
}

/// Trait for node behavior during state computation.
///
/// Implementors can compute new node state based on current state, connected edge flows,
/// external inputs, and time step.
///
pub trait NodeBehavior {
    /// Computes the new state for this node given current conditions.
    ///
    /// # Arguments
    /// * `current_state` - The current runtime state of this node
    /// * `connected_flows` - Flows from all edges connected to this node
    /// * `external_inputs` - External overrides that may apply to this node
    /// * `node_index` - The index of this node in the network
    /// * `dt` - Time step for integration
    ///
    /// # Returns
    /// The updated `CurrentNodeState` for this node
    fn compute_new_state(
        &self,
        current_state: &CurrentNodeState,
        connected_flows: &Vec<(EdgeIndex, FlowResult)>,
        external_inputs: &ExternalInputs,
        node_index: NodeIndex,
        dt: f64,
    ) -> CurrentNodeState;

    /// Returns pressure at a specific connection height for this node type.
    ///
    /// # Arguments
    /// * `state` - Current node state
    /// * `connection_height` - Height from node bottom to connection (m)
    /// * `default_pressure` - Default pressure if node doesn't support height-dependent pressure
    ///
    /// # Returns
    /// Pressure in Pa at the connection height
    fn pressure_at_height(&self, state: &CurrentNodeState, connection_height: f64, default_pressure: f64) -> f64 {
        default_pressure
    }
}

impl NodeBehavior for ThermalNode {
    fn compute_new_state(
        &self,
        current_state: &CurrentNodeState,
        connected_flows: &Vec<(EdgeIndex, FlowResult)>,
        external_inputs: &ExternalInputs,
        node_index: NodeIndex,
        dt: f64,
    ) -> CurrentNodeState {
        // Ensure current state is thermal
        let CurrentNodeState::Thermal {
            temperature: current_temp,
        } = current_state
        else {
            panic!("ThermalNode received non-thermal current state");
        };

        // Calculate net thermal power from connected flows
        let mut net_power = 0.0;
        for (_edge_idx, flow) in connected_flows {
            if let FlowResult::Thermal { thermal_power } = flow {
                // Note: Direction is not properly handled without edge_nodes
                // For now, sum all thermal powers as positive
                net_power += thermal_power;
            }
        }

        // Check for external temperature override
        let new_temp =
            if let Some(&override_temp) = external_inputs.temperature_overrides.get(&node_index) {
                override_temp
            } else {
                // Update temperature: dT = power * dt / mass
                current_temp + net_power * dt / self.thermal_mass
            };

        CurrentNodeState::Thermal {
            temperature: new_temp,
        }
    }
}

impl NodeBehavior for ThermalBoundary {
    fn compute_new_state(
        &self,
        current_state: &CurrentNodeState,
        _connected_flows: &Vec<(EdgeIndex, FlowResult)>,
        _external_inputs: &ExternalInputs,
        _node_index: NodeIndex,
        _dt: f64,
    ) -> CurrentNodeState {
        // Boundary maintains fixed temperature
        // Use current state
        current_state.clone()
    }
}

impl NodeBehavior for FluidBoundary {
    fn compute_new_state(
        &self,
        current_state: &CurrentNodeState,
        _connected_flows: &Vec<(EdgeIndex, FlowResult)>,
        _external_inputs: &ExternalInputs,
        _node_index: NodeIndex,
        _dt: f64,
    ) -> CurrentNodeState {
        // Boundary maintains fixed fluid properties
        // Use current state
        current_state.clone()
    }
}

impl NodeBehavior for FluidJunction {
    fn compute_new_state(
        &self,
        current_state: &CurrentNodeState,
        _connected_flows: &Vec<(EdgeIndex, FlowResult)>,
        _external_inputs: &ExternalInputs,
        _node_index: NodeIndex,
        _dt: f64,
    ) -> CurrentNodeState {
        // Junction is a placeholder, should be replaced during simplification
        // For now, assume it behaves like a constant volume fluid node
        current_state.clone()
    }
}

/// Concrete implementation for thermal edges.

impl NodeBehavior for FluidNode {
    fn compute_new_state(
        &self,
        current_state: &CurrentNodeState,
        connected_flows: &Vec<(EdgeIndex, FlowResult)>,
        external_inputs: &ExternalInputs,
        node_index: NodeIndex,
        dt: f64,
    ) -> CurrentNodeState {
        // Placeholder: Basic fluid node state update
        // For now, just return the current state with minor changes
        // TODO: Implement proper fluid physics (heat transfer, mass balance, etc.)
        let CurrentNodeState::Fluid {
            liquid_temperature: current_liquid_temp,
            gas_temperature: current_gas_temp,
            pressure: current_pressure,
            volume: current_volume,
            liquid_fraction: current_liquid_fraction,
        } = current_state
        else {
            panic!("FluidNode received non-fluid current state");
        };

        // Placeholder: Simple heat accumulation from mass flows
        // This is highly simplified and should be replaced with proper thermodynamics
        let mut liquid_energy_change = 0.0;
        let mut gas_energy_change = 0.0;
        let mut mass_flow_in = 0.0;
        for (_edge_idx, flow) in connected_flows {
            if let FlowResult::Fluid {
                liquid_mass_flow,
                gas_mass_flow,
                liquid_specific_enthalpy,
                gas_specific_enthalpy,
            } = flow
            {
                // Note: Direction is not handled properly without network reference
                liquid_energy_change += liquid_mass_flow * liquid_specific_enthalpy;
                gas_energy_change += gas_mass_flow * gas_specific_enthalpy;
                mass_flow_in += liquid_mass_flow + gas_mass_flow;
            }
        }

        // Placeholder: New temperature as average based on energy
        let total_energy = liquid_energy_change + gas_energy_change;
        let new_liquid_temp =
            current_liquid_temp.unwrap_or(300.0) + total_energy * dt / (1000.0 * current_volume); // Dummy specific heat
        let new_gas_temp =
            current_gas_temp.unwrap_or(300.0) + total_energy * dt / (1000.0 * current_volume);
        let new_pressure = current_pressure + mass_flow_in * dt / current_volume; // Simple pressure from mass
        let new_volume = *current_volume; // Assume constant volume
        let new_liquid_fraction = *current_liquid_fraction; // Placeholder

        CurrentNodeState::Fluid {
            liquid_temperature: Some(new_liquid_temp),
            gas_temperature: Some(new_gas_temp),
            pressure: new_pressure,
            volume: new_volume,
            liquid_fraction: new_liquid_fraction,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_thermal_node_compute_new_state() {
        let node = ThermalNode {
            name: Some("test".to_string()),
            thermal_mass: 100.0,
            initial_temp: 300.0,
        };
        let current_state = CurrentNodeState::Thermal { temperature: 300.0 };
        let connected_flows = vec![(
            EdgeIndex::new(0, 1).unwrap(),
            FlowResult::Thermal { thermal_power: 100.0 },
        )];
        let external_inputs = ExternalInputs {
            temperature_overrides: std::collections::HashMap::new(),
            pressure_overrides: std::collections::HashMap::new(),
        };
        let node_index = NodeIndex::new(0, 1).unwrap();
        let dt = 1.0;

        let new_state = node.compute_new_state(
            &current_state,
            &connected_flows,
            &external_inputs,
            node_index,
            dt,
        );

        match new_state {
            CurrentNodeState::Thermal { temperature } => {
                // Expected: 300 + 100 * 1 / 100 = 301
                assert!((temperature - 301.0).abs() < 1e-6);
            }
            _ => panic!("Wrong state type"),
        }
    }

    #[test]
    fn test_thermal_boundary_no_change() {
        let node = ThermalBoundary {
            name: Some("boundary".to_string()),
            temperature: 200.0,
        };
        let current_state = CurrentNodeState::ThermalBoundary { temperature: 200.0 };
        let connected_flows = vec![];
        let external_inputs = ExternalInputs {
            temperature_overrides: std::collections::HashMap::new(),
            pressure_overrides: std::collections::HashMap::new(),
        };
        let node_index = NodeIndex::new(0, 1).unwrap();
        let dt = 1.0;

        let new_state = node.compute_new_state(
            &current_state,
            &connected_flows,
            &external_inputs,
            node_index,
            dt,
        );

        assert_eq!(new_state, current_state);
    }
}

#[derive(Debug, Clone)]
pub struct ThermalEdgeImpl {
    pub name: Option<String>,
    pub conductance: f64,
    pub emissivity: Option<f64>, // For radiation
}

impl EdgeBehavior for ThermalEdgeImpl {
    fn compute_flow(
        &self,
        node_a_index: NodeIndex,
        node_a_state: &CurrentNodeState,
        connection_height_a: f64,
        node_b_index: NodeIndex,
        node_b_state: &CurrentNodeState,
        connection_height_b: f64,
        external_inputs: &ExternalInputs,
    ) -> FlowResult {
        // Extract temperatures from node states, with external overrides
        let (temp_a, temp_b) = self.get_connected_edge_temps(node_a_state, node_b_state);

        // Calculate thermal power from A to B
        let thermal_power = self.conductance * (temp_a - temp_b);
        FlowResult::Thermal {
            thermal_power: thermal_power,
        }
    }
}

/// Concrete implementation for fluid edges.
#[derive(Debug, Clone)]
pub struct FluidEdgeImpl {
    pub name: Option<String>,
    pub resistance: f64,  // Analogous to conductance
    pub valve_state: f64, // 0..1 for open/close
}

impl EdgeBehavior for FluidEdgeImpl {
    fn compute_flow(
        &self,
        node_a_index: NodeIndex,
        node_a_state: &CurrentNodeState,
        connection_height_a: f64,
        node_b_index: NodeIndex,
        node_b_state: &CurrentNodeState,
        connection_height_b: f64,
        _external_inputs: &ExternalInputs,
    ) -> FlowResult {
        // Get connected fluid states
        let (
            (pressure_a, _temp_a_liquid, _temp_a_gas, _liquid_fraction_a),
            (pressure_b, _temp_b_liquid, _temp_b_gas, _liquid_fraction_b),
        ) = self.get_connected_fluid_states(node_a_state, node_b_state);

        // Placeholder calculation: mass flow based on pressure difference
        let pressure_diff = pressure_a - pressure_b;
        let flow_coeff = 1.0 / (self.resistance + 0.01) * self.valve_state; // Prevent division by zero
        let liquid_mass_flow = flow_coeff * pressure_diff;
        let gas_mass_flow = 0.0; // Placeholder
        FlowResult::Fluid {
            liquid_mass_flow: liquid_mass_flow,
            gas_mass_flow: gas_mass_flow,
            liquid_specific_enthalpy: 4200., // Dummy, should be based on upstream conditions
            gas_specific_enthalpy: 2000.,    // Dummy
        }
    }
}
