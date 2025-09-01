//! Fluid Solver Library
//!
//! This library provides the core implementation for the lumped element thermal/fluid simulation.

pub mod algorithms;
pub mod elements;
pub mod indices;
pub mod network;

// Re-export main types for convenience
pub use algorithms::{LookupTable1D, LookupTable2D};
pub use elements::{Edge, Node};
pub use indices::{EdgeIndex, NodeIndex};
pub use network::FluidNetwork;
