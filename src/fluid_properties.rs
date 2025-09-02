//! Module with fluid properties.
//! We add getters to the properties to hide the underlying implementation.
//! This allows us to increase/decrease fidelity without many adjustments to callsites

use crate::LookupTable1D;
use crate::units::*;

/// Structure with methods to get fluid properties used in this library
pub struct Fluid {
    pub name: String,
    pub molar_mass: f64,
    gas: GasProps,
    liquid: LiqProps,
    specific_gas_const: f64,
    min_pressure: Pressure,
    max_pressure: Pressure,
    min_temperature: Temperature,
    max_temperature: Temperature,
    pressure_to_sat_temp_table: LookupTable1D,
    pressure_to_heat_of_vaporization_table: LookupTable1D,
}

struct GasProps {
    ref_temp: Temperature,
    ref_pressure: Pressure,
    ref_density: Density,
    ref_cp: f64,
    ref_cv: f64,
    ref_specific_enthalpy: f64,
    ref_compressibility: f64,
}

struct LiqProps {
    ref_temp: Temperature,
    ref_pressure: Pressure,
    ref_density: Density,
    ref_cp: f64,
    ref_cv: f64,
    ref_specific_enthalpy: f64,
}

impl Fluid {
    /// Get the specific gas constant for this fluid
    pub fn specific_gas_const(&self) -> f64 {
        self.specific_gas_const
    }
    /// Get the saturation temperature at a given pressure
    pub fn saturation_temperature(&self, press: Pressure) -> Temperature {
        self.pressure_to_sat_temp_table
            .interpolate_linear_stateful(press)
    }
    /// Get the heat fo vaporization at a given pressure
    pub fn heat_of_vaporization(&self, press: Pressure) -> SpecificEnergy {
        self.pressure_to_heat_of_vaporization_table
            .interpolate_linear_stateful(press)
    }
    /// Get pressure given mass, volume and temperature
    pub fn calc_pressure(&self, mass: Mass, vol: Volume, temp: Temperature) -> Pressure {
        // For now use calorically perfect gas with non unity compressibility
        (self.gas.ref_compressibility * mass * self.specific_gas_const * temp / vol)
            .clamp(self.min_pressure, self.max_pressure)
    }
    /// Get temperature from state
    pub fn calc_temperature(
        &self,
        mass: Mass,
        energy: Energy,
        pressure: Pressure,
        vol: Volume,
    ) -> Temperature {
    }

    pub fn temp_from_enthalpy(&self, enthalpy: f64) -> Temperature {
        (enthalpy - self.gas.ref_specific_enthalpy) / self.gas.ref_cp + self.gas.ref_temp
    }

    pub fn enthalpy_from_temp(&self, temp: Temperature) -> f64 {
        (temp - self.gas.ref_temp) * self.gas.ref_cp + self.gas.ref_specific_enthalpy
    }
}
