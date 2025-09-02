// Base units for US customary system
// Type aliases for fundamental units using f64
pub type Time = f64; // seconds
pub type Length = f64; // feet
pub type Mass = f64; // slug
pub type Force = f64; // pound-force (lbf)
pub type Temperature = f64; // Rankine

// Derived unit types (not base, but commonly used)
pub type Pressure = f64; // lbf/ft² (base pressure unit)
pub type Area = f64; // ft²
pub type Volume = f64; // ft³
pub type MassFlow = f64; // slug/s
pub type Energy = f64; // ft*lbf
pub type SpecificEnergy = f64; // ft*lbf/slug
pub type Power = f64; // ft*lbf/s
pub type Density = f64; // ft*lbf/s

// Common physical constants
pub const GRAVITY_ACCELERATION: f64 = 32.174; // ft/s² (standard gravitational acceleration)
pub const WATER_DENSITY: f64 = 1.934593; // slug/ft³ (density of water)

// Standard atmospheric conditions
pub const ATMOSPHERIC_PRESSURE_PSI: f64 = 14.696; // psi (standard atmospheric pressure)
pub const ATMOSPHERIC_PRESSURE_BASE: Pressure = ATMOSPHERIC_PRESSURE_PSI * 144.0; // lbf/ft²
pub const STANDARD_AMBIENT_TEMPERATURE_F: f64 = 59.0; // °F (standard temperature for atmospheric pressure)
pub const STANDARD_AMBIENT_TEMPERATURE_BASE: Temperature = STANDARD_AMBIENT_TEMPERATURE_F + 459.67; // Rankine

// Conversion functions from common units to base units

/// Convert pressure from psi (pounds per square inch) to lbf/ft²
///
/// # Arguments
/// * `psi` - Pressure in psi
///
/// # Returns
/// Pressure in lbf/ft²
///
/// # Examples
/// ```
/// let base_pressure = psi_to_psf(14.7); // Atmospheric pressure
/// ```
pub fn psi_to_psf(psi: f64) -> Pressure {
    psi * 144.0
}

/// Convert pressure from lbf/ft² to psi
///
/// # Arguments
/// * `psf` - Pressure in lbf/ft²
///
/// # Returns
/// Pressure in psi
pub fn psf_to_psi(psf: Pressure) -> f64 {
    psf / 144.0
}

// Temperature conversion functions (Rankine is the base unit)
/// Convert Fahrenheit to Rankine
pub fn fahrenheit_to_rankine(f: f64) -> Temperature {
    f + 459.67
}

/// Convert Rankine to Fahrenheit
pub fn rankine_to_fahrenheit(r: Temperature) -> f64 {
    r - 459.67
}

/// Convert Celsius to Rankine
pub fn celsius_to_rankine(c: f64) -> Temperature {
    (c + 273.15) * 9.0 / 5.0
}

/// Convert Rankine to Celsius
pub fn rankine_to_celsius(r: Temperature) -> f64 {
    (r * 5.0 / 9.0) - 273.15
}

/// Convert Kelvin to Rankine
pub fn kelvin_to_rankine(k: f64) -> Temperature {
    k * 9.0 / 5.0
}

/// Convert Rankine to Kelvin
pub fn rankine_to_kelvin(r: Temperature) -> f64 {
    r * 5.0 / 9.0
}
