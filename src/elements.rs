//! Defines the core elements of the fluid network.
//!
//! This module contains structs and enums for nodes, edges, and their properties.

use crate::indices::NodeIndex;
use std::collections::HashMap;

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
}
/// Represents external inputs that can override node states at runtime.
#[derive(Debug, Clone)]
pub struct ExternalInputs {
    pub temperature_overrides: HashMap<NodeIndex, f64>,
    pub pressure_overrides: HashMap<NodeIndex, f64>,
}

/// Represents the current runtime state of a node.
/// Separate temperatures for liquid and gas phases in fluid nodes.
/// Sign conventions: Flows are calculated from node A to node B, positive means from A to B.
#[derive(Debug, Clone)]
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


/// Trait for edge behavior, including computing flow.
pub trait EdgeBehavior {
    fn compute_flow(
        &self,
        node_a_index: NodeIndex,
        node_a_state: &CurrentNodeState,
        node_b_index: NodeIndex,
        node_b_state: &CurrentNodeState,
        external_inputs: &ExternalInputs,
    ) -> FlowResult;
}

/// Concrete implementation for thermal edges.
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
        node_b_index: NodeIndex,
        node_b_state: &CurrentNodeState,
        external_inputs: &ExternalInputs,
    ) -> FlowResult {
        // Extract temperatures from node states, with external overrides
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
        node_b_index: NodeIndex,
        node_b_state: &CurrentNodeState,
        external_inputs: &ExternalInputs,
    ) -> FlowResult {
        let (mut pressure_a, mut temp_a_liquid, mut temp_a_gas) = match node_a_state {
            CurrentNodeState::Fluid {
                liquid_temperature,
                gas_temperature,
                pressure,
                ..
            } => (*pressure, *liquid_temperature, *gas_temperature),
            CurrentNodeState::FluidBoundary {
                liquid_temperature,
                gas_temperature,
                pressure,
                liquid_fraction,
            } => (*pressure, *liquid_temperature, *gas_temperature),
            _ => {
                unreachable!("Validation should ensure all edge connecitons are valid.");
            } // Invalid for fluid edges
        };
        let (mut pressure_b, mut temp_b_liquid, mut temp_b_gas) = match node_b_state {
            CurrentNodeState::Fluid {
                liquid_temperature,
                gas_temperature,
                pressure,
                ..
            } => (*pressure, *liquid_temperature, *gas_temperature),
            CurrentNodeState::FluidBoundary {
                liquid_temperature,
                gas_temperature,
                pressure,
                liquid_fraction,
            } => (*pressure, *liquid_temperature, *gas_temperature),
            _ => {
                unreachable!("Validation should ensure all edge connecitons are valid.");
            } // Invalid for fluid edges
        };

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

/// Refactored Edge struct that uses trait objects for flexibility.
pub struct Edge {
    pub name: Option<String>,
    pub behavior: Box<dyn EdgeBehavior>,
    pub edge_type: EdgeType,
}

impl std::fmt::Debug for Edge {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Edge {{ ... }}")
    }
}
