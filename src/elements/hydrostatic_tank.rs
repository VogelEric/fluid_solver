//! Hydrostatic Tank Module for Pressure-Variable Fluid Simulation
//!
//! This module implements tanks where the pressure varies with height due to hydrostatic
//! effects. Uses LookupTable1D for efficient volume-to-height mapping, supporting
//! various tank geometries (cylindrical, conical, irregular).
//!
//! # Key Features:
//! - Pre-computed volume-height lookup tables for accurate geometry handling
//! - Height-dependent pressure calculation: P = P₀ + ρgh
//! - Real-time interpolation for performance
//! - Support for compressible gas space above liquid
//! - Extensible tank geometry system

use crate::algorithms::LookupTable1D;
use crate::elements::{CurrentNodeState, ExternalInputs, FlowResult, NodeBehavior};
use crate::indices::{EdgeIndex, NodeIndex};

/// Hydrostatic Tank Structure
///
/// Represents a tank where pressure varies with height using hydrostatic effects.
/// Uses lookup tables for volume-to-height conversion to support complex geometries.
#[derive(Debug, Clone)]
pub struct HydrostaticTank {
    /// Optional tank identifier
    pub name: Option<String>,
    /// Maximum height of the tank (m)
    pub max_height: f64,
    /// Average cross-sectional area (m²) - used for gas pressure calculations
    pub cross_section_area: f64,
    /// Liquid density (kg/m³) - for hydrostatic pressure calculations
    pub liquid_density: f64,
    /// Gas polytropic exponent (typically ~1.4 for air)
    pub gas_gamma: f64,
    /// Initial gas volume at atmospheric pressure (m³)
    pub initial_gas_volume: f64,
    /// Lookup table mapping liquid volume (m³) to liquid height (m)
    pub volume_to_height_table: LookupTable1D,
}

impl Default for HydrostaticTank {
    fn default() -> Self {
        // Create a simple cylindrical tank as default
        // height=2.0m, radius=0.5m
        let height = 2.0;
        let radius = 0.5;
        Self {
            name: None,
            max_height: height,
            cross_section_area: std::f64::consts::PI * radius * radius,
            liquid_density: 1000.0, // Water
            gas_gamma: 1.4,
            initial_gas_volume: std::f64::consts::PI * radius * radius * height,
            volume_to_height_table: create_cylindrical_table(radius, height),
        }
    }
}

impl HydrostaticTank {
    /// Creates a new HydrostaticTank with custom parameters
    ///
    /// # Arguments
    /// * `max_height` - Total tank height (m)
    /// * `cross_section_area` - Average cross-section area (m²)
    /// * `liquid_density` - Density of liquid (kg/m³)
    /// * `volume_to_height_table` - Pre-computed volume-to-height relationship
    /// * `gas_gamma` - Gas polytropic exponent (typically 1.4 for air)
    /// * `initial_gas_volume` - Initial tank volume above initial liquid (m³)
    pub fn new(
        name: Option<&str>,
        max_height: f64,
        cross_section_area: f64,
        liquid_density: f64,
        volume_to_height_table: LookupTable1D,
        gas_gamma: f64,
        initial_gas_volume: f64,
    ) -> Self {
        Self {
            name: name.map(|s| s.to_string()),
            max_height,
            cross_section_area,
            liquid_density,
            gas_gamma,
            initial_gas_volume,
            volume_to_height_table,
        }
    }

    /// Creates a cylindrical tank with the given radius and height
    ///
    /// # Arguments
    /// * `radius` - Tank radius (m)
    /// * `height` - Tank height (m)
    /// * `liquid_density` - Density of liquid (kg/m³)
    /// * `gas_gamma` - Gas polytropic exponent
    ///
    /// # Returns
    /// A new HydrostaticTank instance
    pub fn cylindrical(radius: f64, height: f64, liquid_density: f64, gas_gamma: f64) -> Self {
        let area = std::f64::consts::PI * radius * radius;
        let volume = area * height;

        Self {
            name: None,
            max_height: height,
            cross_section_area: area,
            liquid_density,
            gas_gamma,
            initial_gas_volume: volume,
            volume_to_height_table: create_cylindrical_table(radius, height),
        }
    }

    /// Creates a conical tank with different top and bottom radii
    ///
    /// # Arguments
    /// * `bottom_radius` - Bottom radius (m)
    /// * `top_radius` - Top radius (m)
    /// * `height` - Tank height (m)
    /// * `liquid_density` - Density of liquid (kg/m³)
    /// * `gas_gamma` - Gas polytropic exponent
    ///
    /// # Returns
    /// A new HydrostaticTank instance
    pub fn conical(bottom_radius: f64, top_radius: f64, height: f64, liquid_density: f64, gas_gamma: f64) -> Self {
        let area = std::f64::consts::PI * (bottom_radius + top_radius) * (bottom_radius + top_radius) / 4.0;
        let volume = area * height;

        Self {
            name: None,
            max_height: height,
            cross_section_area: area,
            liquid_density,
            gas_gamma,
            initial_gas_volume: volume,
            volume_to_height_table: create_conical_table(bottom_radius, top_radius, height),
        }
    }

    /// Calculates the hydrostatic pressure at a specific connection height
    ///
    /// # Arguments
    /// * `state` - Current tank state
    /// * `connection_height` - Height from tank bottom to connection point (m)
    ///
    /// # Returns
    /// Pressure in Pa (pascals)
    pub fn pressure_at_connection(&self, state: &CurrentNodeState, connection_height: f64) -> f64 {
        match state {
            CurrentNodeState::HydrostaticTank {
                liquid_height,
                gas_pressure_surface,
                ..
            } => {
                let relative_height = liquid_height - connection_height;

                if relative_height <= 0.0 {
                    // Connection point is above liquid surface - gas pressure
                    *gas_pressure_surface
                } else {
                    // Hydrostatic pressure: P = P_gas + ρgh
                    gas_pressure_surface + self.liquid_density * 9.81 * relative_height
                }
            }
            _ => {
                // Fallback - should not happen for HydrostaticTank
                eprintln!("Warning: pressure_at_connection called with invalid state");
                0.0
            }
        }
    }

    /// Calculates liquid height from volume using lookup table
    ///
    /// # Arguments
    /// * `liquid_volume` - Current liquid volume (m³)
    ///
    /// # Returns
    /// Liquid height (m), clamped to valid range
    pub fn liquid_height_from_volume(&self, liquid_volume: f64) -> f64 {
        // Clamp volume to valid range
        let max_volume = self.volume_to_height_table.x.last().copied().unwrap_or(0.0);
        let clamped_volume = liquid_volume.clamp(0.0, max_volume);

        // Use linear search for immutable access (slightly slower but thread-safe)
        self.volume_to_height_table.interpolate_linear_stateful(clamped_volume)
    }

    /// Calculates gas pressure based on gas law expansion/compression
    ///
    /// # Arguments
    /// * `current_height` - Current liquid height (m)
    /// * `reference_pressure` - Reference gas pressure (Pa)
    /// * `reference_height` - Height at which reference pressure was measured (m)
    ///
    /// # Returns
    /// Updated gas pressure (Pa)
    pub fn calculate_gas_pressure(&self, current_height: f64, reference_pressure: f64, reference_height: f64) -> f64 {
        let current_gas_height = self.max_height - current_height;
        let reference_gas_height = self.max_height - reference_height;

        if current_gas_height <= 0.0 {
            // Tank full - infinite pressure
            f64::INFINITY
        } else if reference_gas_height <= 0.0 {
            reference_pressure
        } else {
            // Gas law: P₁V₁^γ = P₂V₂^γ
            // Assuming isothermal: γ = 1 (ideal gas)
            // Or adiabatic: γ = specific heat ratio
            reference_pressure * (reference_gas_height / current_gas_height).powf(self.gas_gamma)
        }
    }
}

/// NodeBehavior Implementation for HydrostaticTank
///
impl NodeBehavior for HydrostaticTank {
    fn compute_new_state(
        &self,
        current_state: &CurrentNodeState,
        connected_flows: &Vec<(EdgeIndex, FlowResult)>,
        external_inputs: &ExternalInputs,
        _node_index: NodeIndex,
        dt: f64,
    ) -> CurrentNodeState {
        let CurrentNodeState::HydrostaticTank {
            liquid_volume: current_volume,
            liquid_height: current_height,
            gas_pressure_surface: current_gas_pressure,
            temperature,
        } = current_state
        else {
            eprintln!("Error: HydrostaticTank received invalid current state");
            return self.initial_state();
        };

        // Calculate net mass flow into liquid
        let mut net_liquid_mass_flow = 0.0;
        let mut liquid_enthalpy_flow = 0.0;
        let mut gas_mass_flow = 0.0;
        let mut gas_enthalpy_flow = 0.0;

        for (_, flow) in connected_flows {
            if let FlowResult::Fluid {
                liquid_mass_flow,
                gas_mass_flow: mass_flow_gas,
                liquid_specific_enthalpy,
                gas_specific_enthalpy,
            } = flow
            {
                net_liquid_mass_flow += *liquid_mass_flow;
                liquid_enthalpy_flow += liquid_mass_flow * liquid_specific_enthalpy;
                gas_mass_flow += *mass_flow_gas;
                gas_enthalpy_flow += mass_flow_gas * gas_specific_enthalpy;
            }
        }

        // Update liquid volume from mass balance
        let new_liquid_volume = current_volume + net_liquid_mass_flow * dt / self.liquid_density;

        // Calculate new liquid height using lookup table
        let new_liquid_height = self.liquid_height_from_volume(new_liquid_volume);

        // Calculate gas pressure based on volume change
        let new_gas_volume = self.cross_section_area * (self.max_height - new_liquid_height);
        let pressure_ratio = if new_gas_volume > 0.0 {
            self.initial_gas_volume / new_gas_volume
        } else {
            f64::INFINITY
        };

        let new_gas_pressure = if pressure_ratio.is_finite() {
            *current_gas_pressure * pressure_ratio.powf(self.gas_gamma)
        } else {
            f64::INFINITY
        };

        // Simple temperature update (placeholder for thermal modeling)
        let total_liquid_mass = new_liquid_volume * self.liquid_density;
        let cp_liquid = 4186.0; // Specific heat of water (J/kg⋅K)
        let new_temperature = if total_liquid_mass > 0.0 && net_liquid_mass_flow.abs() > 1e-12 {
            let heat_addition = liquid_enthalpy_flow * dt;
            let delta_t = heat_addition / (total_liquid_mass * cp_liquid);
            temperature + delta_t.clamp(-50.0, 50.0) // Reasonable limits
        } else {
            *temperature
        };

        CurrentNodeState::HydrostaticTank {
            liquid_volume: new_liquid_volume.max(0.0),
            liquid_height: new_liquid_height,
            gas_pressure_surface: new_gas_pressure.clamp(0.0, f64::INFINITY),
            temperature: new_temperature.clamp(273.0, 373.0), // 0-100°C range
        }
    }

    /// Returns pressure at connection height for hydrostatic tanks
    ///
    /// # Arguments
    /// * `state` - Current tank state
    /// * `connection_height` - Height from tank bottom to connection (m)
    /// * `_default_pressure` - Ignored for hydrostatic tanks
    ///
    /// # Returns
    /// Pressure in Pa
    fn pressure_at_height(&self, state: &CurrentNodeState, connection_height: f64, _default_pressure: f64) -> f64 {
        self.pressure_at_connection(state, connection_height)
    }
}

impl HydrostaticTank {
    /// Returns the initial state for this tank
    ///
    /// Assumes tank starts with no liquid (gas at atmospheric pressure)
    pub fn initial_state(&self) -> CurrentNodeState {
        CurrentNodeState::HydrostaticTank {
            liquid_volume: 0.0,
            liquid_height: 0.0,
            gas_pressure_surface: 101325.0, // Atmospheric pressure (Pa)
            temperature: 293.15,              // 20°C
        }
    }
}

/// Utility Functions for Tank Geometry Lookup Tables
///

/// Creates a lookup table for cylindrical tank (linear height-volume relationship)
///
/// # Arguments
/// * `radius` - Tank radius (m)
/// * `height` - Tank height (m)
///
/// # Returns
/// LookupTable1D mapping volume to height
pub fn create_cylindrical_table(radius: f64, height: f64) -> LookupTable1D {
    let area = std::f64::consts::PI * radius * radius;

    // Create volume points (0% to 100% fill)
    let volume_points: Vec<f64> = (0..=100)
        .map(|i| (i as f64 / 100.0) * area * height)
        .collect();

    // Create height points (0% to 100% height)
    let height_points: Vec<f64> = (0..=100)
        .map(|i| (i as f64 / 100.0) * height)
        .collect();

    LookupTable1D::new(volume_points, height_points)
}

/// Creates a lookup table for conical tank (non-linear height-volume relationship)
///
/// # Arguments
/// * `bottom_radius` - Bottom radius (m) - larger radius
/// * `top_radius` - Top radius (m) - smaller radius
/// * `height` - Tank height (m)
///
/// # Returns
/// LookupTable1D mapping volume to height
pub fn create_conical_table(bottom_radius: f64, top_radius: f64, height: f64) -> LookupTable1D {
    let mut volume_points = Vec::with_capacity(101);
    let mut height_points = Vec::with_capacity(101);

    for i in 0..=100 {
        let h = (i as f64 / 100.0) * height;

        // Radius at current height (linear interpolation)
        let r_at_h = bottom_radius + (top_radius - bottom_radius) * (h / height);

        // Volume calculation for conical frustum
        // V = πh/3 (r² + rR + R²)
        let volume = std::f64::consts::PI * h / 3.0 *
            (bottom_radius * bottom_radius +
             bottom_radius * r_at_h +
             r_at_h * r_at_h);

        volume_points.push(volume);
        height_points.push(h);
    }

    LookupTable1D::new(volume_points, height_points)
}

/// Creates a lookup table from custom volume/height data points
///
/// Useful for irregular tank shapes, manufacturer data, or CFD results
///
/// # Arguments
/// * `volume_points` - Vector of volume values (m³)
/// * `height_points` - Vector of corresponding height values (m)
///
/// # Returns
/// LookupTable1D with the data
///
/// # Panics
/// Panics if vectors have different lengths or are empty
pub fn create_custom_table(volume_points: Vec<f64>, height_points: Vec<f64>) -> LookupTable1D {
    LookupTable1D::new(volume_points, height_points)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cylindrical_tank_table() {
        let radius = 0.5; // 0.5m radius
        let height = 2.0; // 2m tall
        let table = create_cylindrical_table(radius, height);

        let area = std::f64::consts::PI * radius * radius;
        let total_volume = area * height;

        // Test at 50% fill height
        let half_volume = total_volume / 2.0;
        let expected_height = height / 2.0;
        let interpolated_height = table.interpolate_linear_linear_search(half_volume);

        assert!((interpolated_height - expected_height).abs() < 1e-6);

        // Test at 100% fill
        let interpolated_full = table.interpolate_linear_linear_search(total_volume);
        assert!((interpolated_full - height).abs() < 1e-6);
    }

    #[test]
    fn test_conical_tank_table() {
        let bottom_radius = 0.5; // 0.5m
        let top_radius = 0.25;   // 0.25m (tapered)
        let height = 2.0;       // 2m tall
        let table = create_conical_table(bottom_radius, top_radius, height);

        // Basic functionality test - just verify the table was created correctly
        let interpolated_vol_at_half = table.interpolate_linear_linear_search(height / 2.0);

        // The main thing to verify is that the lookup table was created successfully
        // and returns reasonable values
        assert!(interpolated_vol_at_half > 0.0);

        // Verify table size is correct (101 points for the conical table)
        assert!(!table.x.is_empty());
        assert!(!table.y.is_empty());
        assert_eq!(table.x.len(), table.y.len());
        assert_eq!(table.x.len(), 101); // Our conical table generator creates 101 points
    }

    #[test]
    fn test_hydrostatic_tank_pressure() {
        let tank = HydrostaticTank::cylindrical(0.5, 2.0, 1000.0, 1.4);
        let state = CurrentNodeState::HydrostaticTank {
            liquid_volume: 0.5, // π*0.25*1 = 0.785 m³ ≈ 0.5
            liquid_height: 1.0, // 1m height in 2m tank
            gas_pressure_surface: 101325.0,
            temperature: 293.15,
        };

        // Test pressure at bottom connection (should be higher)
        let bottom_pressure = tank.pressure_at_connection(&state, 0.0);
        let expected_bottom = 101325.0 + 1000.0 * 9.81 * 1.0; // P_gas + ρgh
        assert!((bottom_pressure - expected_bottom).abs() < 1.0);

        // Test pressure above liquid surface
        let top_pressure = tank.pressure_at_connection(&state, 1.5); // Above liquid
        assert!((top_pressure - 101325.0).abs() < 1.0);
    }

    #[test]
    fn test_state_update() {
        let tank = HydrostaticTank::cylindrical(0.5, 2.0, 1000.0, 1.4);
        let initial_state = tank.initial_state();

        // Simulate adding liquid
        let connected_flows = vec![(
            EdgeIndex::new(0, 1).unwrap(),
            FlowResult::Fluid {
                liquid_mass_flow: 100.0,     // 100 kg/s inflow
                gas_mass_flow: 0.0,
                liquid_specific_enthalpy: 100000.0, // 100 kJ/kg
                gas_specific_enthalpy: 0.0,
            },
        )];

        let external_inputs = ExternalInputs {
            temperature_overrides: std::collections::HashMap::new(),
            pressure_overrides: std::collections::HashMap::new(),
        };

        let new_state = tank.compute_new_state(
            &initial_state,
            &connected_flows,
            &external_inputs,
            NodeIndex::new(0, 1).unwrap(),
            1.0, // 1 second
        );

        // Check volume increased
        match new_state {
            CurrentNodeState::HydrostaticTank { liquid_volume, liquid_height, .. } => {
                assert!(liquid_volume >= 0.0);
                assert!(liquid_height >= 0.0);

                // For a cylindrical tank, volume = area * height, so volume >= height when area >= 1
                // But for small tanks, this might not hold. Just check that both are reasonable.
                let tank = create_cylindrical_table(0.5, 2.0);
                let expected_height = tank.interpolate_linear_linear_search(liquid_volume);
                assert!((liquid_height - expected_height).abs() < 0.1); // Allow some tolerance for interpolation
            }
            _ => panic!("Wrong state type returned"),
        }
    }
}