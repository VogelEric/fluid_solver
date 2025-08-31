//! Defines the core elements of the fluid network.
//!
//! This module contains structs and enums for nodes, edges, and their properties.

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

/// Result of edge flow computation.
#[derive(Debug, Clone)]
pub struct FlowResult {
    pub thermal_power: Option<f64>,
    pub liquid_mass_flow: Option<f64>,
    pub gas_mass_flow: Option<f64>,
    pub liquid_specific_enthalpy: Option<f64>,
    pub gas_specific_enthalpy: Option<f64>,
}

/// Trait for edge behavior, including computing flow.
pub trait EdgeBehavior {
    fn compute_flow(&self) -> FlowResult;
}

/// Concrete implementation for thermal edges.
#[derive(Debug, Clone)]
pub struct ThermalEdgeImpl {
    pub name: Option<String>,
    pub conductance: f64,
    pub emissivity: Option<f64>, // For radiation
}

impl EdgeBehavior for ThermalEdgeImpl {
    fn compute_flow(&self) -> FlowResult {
        // Placeholder: in real implementation, calculate power based on temperature difference
        FlowResult {
            thermal_power: Some(self.conductance * 100.0), // Dummy calculation
            liquid_mass_flow: None,
            gas_mass_flow: None,
            liquid_specific_enthalpy: None,
            gas_specific_enthalpy: None,
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
    fn compute_flow(&self) -> FlowResult {
        // Placeholder: in real implementation, calculate mass flows and enthalpies
        FlowResult {
            thermal_power: None,
            liquid_mass_flow: Some(self.resistance * self.valve_state * 10.0), // Dummy
            gas_mass_flow: Some(0.0),
            liquid_specific_enthalpy: Some(4200.0), // Dummy for water
            gas_specific_enthalpy: Some(2000.0),
        }
    }
}

/// Refactored Edge struct that uses trait objects for flexibility.
pub struct Edge {
    pub name: Option<String>,
    pub behavior: Box<dyn EdgeBehavior>,
}

impl std::fmt::Debug for Edge {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Edge {{ ... }}")
    }
}
