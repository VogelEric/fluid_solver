//! Provides custom index types for safe node and edge references.
//!
//! These newtypes prevent accidental index errors and enable validation.

/// Represents a valid index into the nodes vector.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NodeIndex(usize);

impl NodeIndex {
    /// Creates a new NodeIndex if the index is within bounds.
    pub fn new(index: usize, node_count: usize) -> Option<Self> {
        if index < node_count {
            Some(Self(index))
        } else {
            None
        }
    }

    /// Returns the underlying usize value.
    pub fn index(&self) -> usize {
        self.0
    }
}

/// Represents a valid index into the edges vector.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EdgeIndex(usize);

impl EdgeIndex {
    /// Creates a new EdgeIndex if the index is within bounds.
    pub fn new(index: usize, edge_count: usize) -> Option<Self> {
        if index < edge_count {
            Some(Self(index))
        } else {
            None
        }
    }

    /// Returns the underlying usize value.
    pub fn index(&self) -> usize {
        self.0
    }
}
