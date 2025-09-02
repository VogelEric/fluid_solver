//! Algorithms module for common utilities like lookup tables.
//!
//! This module provides 1D and 2D lookup tables with linear interpolation
//! and different search algorithms for breakpoint finding.

use std::cell::Cell;
use std::cmp::Ordering;

/// 1D Lookup Table for fast interpolation of y values from x inputs.
/// Assumes x is sorted in ascending order.
#[derive(Debug, Clone)]
pub struct LookupTable1D {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    last_index: Cell<usize>, // For stateful search optimization
}

impl LookupTable1D {
    /// Creates a new LookupTable1D from x and y vectors.
    /// Sorts the x values if not already sorted.
    ///
    /// # Panics
    /// Panics if x and y have different lengths or x is empty.
    pub fn new(mut x: Vec<f64>, mut y: Vec<f64>) -> Self {
        assert_eq!(x.len(), y.len(), "x and y must have the same length");
        assert!(!x.is_empty(), "x cannot be empty");

        // Sort x and y together
        let mut x_sorted: Vec<f64> = x.iter().cloned().collect();
        let mut y_sorted: Vec<f64> = y.iter().cloned().collect();
        let mut indices: Vec<usize> = (0..x.len()).collect();

        indices.sort_by(|&i, &j| x[i].partial_cmp(&x[j]).unwrap());

        for (new_idx, &old_idx) in indices.iter().enumerate() {
            x_sorted[new_idx] = x[old_idx];
            y_sorted[new_idx] = y[old_idx];
        }

        x = x_sorted;
        y = y_sorted;

        LookupTable1D {
            x,
            y,
            last_index: Cell::new(0),
        }
    }

    /// Linearly interpolates y value for given x using linear search (stateless).
    /// Uses clamping for out-of-bounds values.
    pub fn interpolate_linear_linear_search(&self, input: f64) -> f64 {
        let len = self.x.len();
        if len == 1 {
            return self.y[0];
        }

        if input <= self.x[0] {
            return self.y[0];
        } else if input >= self.x[len - 1] {
            return self.y[len - 1];
        }

        // Linear search for bracket
        for i in 0..len - 1 {
            if self.x[i] <= input && input <= self.x[i + 1] {
                let frac = (input - self.x[i]) / (self.x[i + 1] - self.x[i]);
                return self.y[i] + frac * (self.y[i + 1] - self.y[i]);
            }
        }
        unreachable!(); // Should not reach here due to clamping
    }

    /// Linearly interpolates y value for given x using binary search (stateless).
    /// Uses clamping for out-of-bounds values.
    pub fn interpolate_linear_binary_search(&self, input: f64) -> f64 {
        let len = self.x.len();
        if len == 1 {
            return self.y[0];
        }

        if input <= self.x[0] {
            return self.y[0];
        } else if input >= self.x[len - 1] {
            return self.y[len - 1];
        }

        // Binary search for bracket using standard library
        match self
            .x
            .binary_search_by(|&val| val.partial_cmp(&input).unwrap_or(Ordering::Equal))
        {
            Ok(idx) => self.y[idx],
            Err(insert_idx) => {
                let i = insert_idx.saturating_sub(1);
                let j = i + 1;
                let frac = (input - self.x[i]) / (self.x[j] - self.x[i]);
                self.y[i] + frac * (self.y[j] - self.y[i])
            }
        }
    }

    /// Linearly interpolates y value for given x using stateful search.
    /// Remembers the last index to optimize successive lookups.
    /// Uses clamping for out-of-bounds values.
    pub fn interpolate_linear_stateful(&self, input: f64) -> f64 {
        let len = self.x.len();
        if len == 1 {
            return self.y[0];
        }

        if input <= self.x[0] {
            self.last_index.set(0);
            return self.y[0];
        } else if input >= self.x[len - 1] {
            self.last_index.set(len - 1);
            return self.y[len - 1];
        }

        // Check if near last index (within a few indices)
        let mut i = self.last_index.get().saturating_sub(2);
        let mut found = false;
        for k in 0..5 {
            // Check up to 5 indices around last_index
            let idx = (i + k).min(len - 1);
            if idx >= len - 1 {
                break;
            }
            if self.x[idx] <= input && input <= self.x[idx + 1] {
                i = idx;
                found = true;
                self.last_index.set(i);
                break;
            }
        }

        if !found {
            // Fall back to binary search
            match self
                .x
                .binary_search_by(|&val| val.partial_cmp(&input).unwrap_or(Ordering::Equal))
            {
                Ok(idx) => {
                    self.last_index.set(idx);
                    self.y[idx]
                }
                Err(insert_idx) => {
                    i = insert_idx.saturating_sub(1);
                    self.last_index.set(i);
                    let j = i + 1;
                    let frac = (input - self.x[i]) / (self.x[j] - self.x[i]);
                    self.y[i] + frac * (self.y[j] - self.y[i])
                }
            }
        } else {
            let j = (i + 1).min(len - 1);
            let frac = (input - self.x[i]) / (self.x[j] - self.x[i]);
            self.y[i] + frac * (self.y[j] - self.y[i])
        }
    }
}

/// 2D Lookup Table for bilinear interpolation.
/// Assumes x and y are sorted in ascending order.
#[derive(Debug, Clone)]
pub struct LookupTable2D {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub values: Vec<Vec<f64>>, // values[y_index][x_index]
    last_x_index: Cell<usize>,
    last_y_index: Cell<usize>,
}

impl LookupTable2D {
    /// Creates a new LookupTable2D from x, y, and values vectors.
    /// Sorts x and y, and transposes values accordingly.
    ///
    /// # Panics
    /// Panics if dimensions don't match or are empty.
    #[allow(unused_mut)] // We may decide to mutably use this in the future
    pub fn new(mut x: Vec<f64>, mut y: Vec<f64>, mut values: Vec<Vec<f64>>) -> Self {
        assert!(!x.is_empty() && !y.is_empty(), "x and y cannot be empty");
        assert_eq!(values.len(), y.len(), "values outer length must match y");
        for row in &values {
            assert_eq!(row.len(), x.len(), "values inner lengths must match x");
        }

        // Sort x and y, rearrange values
        let mut x_pairs: Vec<_> = x.iter().cloned().enumerate().collect();
        x_pairs.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        x = x_pairs.iter().map(|(_, val)| *val).collect();

        let mut y_pairs: Vec<_> = y.iter().cloned().enumerate().collect();
        y_pairs.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        y = y_pairs.iter().map(|(_, val)| *val).collect();

        // Rearrange values[y_old][x_old] to values[y_new][x_new]
        let mut new_values = vec![vec![0.0; x.len()]; y.len()];
        for (new_y, (old_y_idx, _)) in y_pairs.iter().enumerate() {
            let old_row = &values[*old_y_idx];
            for (new_x, (old_x_idx, _)) in x_pairs.iter().enumerate() {
                new_values[new_y][new_x] = old_row[*old_x_idx];
            }
        }

        LookupTable2D {
            x,
            y,
            values: new_values,
            last_x_index: Cell::new(0),
            last_y_index: Cell::new(0),
        }
    }

    /// Performs bilinear interpolation using stateful accelerated search.
    /// Uses clamping for out-of-bounds values and remembers last indices.
    pub fn interpolate_bilinear_stateful(&self, x_input: f64, y_input: f64) -> f64 {
        let x_len = self.x.len();
        let y_len = self.y.len();

        if x_len == 1 && y_len == 1 {
            return self.values[0][0];
        }

        // Find x indices
        let x_idx = if x_input <= self.x[0] {
            self.last_x_index.set(0);
            0
        } else if x_input >= self.x[x_len - 1] {
            self.last_x_index.set(x_len - 1);
            x_len - 1
        } else {
            // Try accelerated search around last_x_index
            let start_x = self.last_x_index.get().saturating_sub(2);
            let mut found_x = false;
            let mut x_i = start_x;
            let max_x = x_len - 1;
            for k in 0..5 {
                let idx = (start_x + k).min(max_x);
                if idx == max_x {
                    break;
                }
                if self.x[idx] <= x_input && x_input <= self.x[idx + 1] {
                    x_i = idx;
                    found_x = true;
                    break;
                }
            }
            if found_x {
                x_i
            } else {
                self.find_bracket_binary(&self.x, x_input)
            }
        };

        // Find y indices
        let y_idx = if y_input <= self.y[0] {
            self.last_y_index.set(0);
            0
        } else if y_input >= self.y[y_len - 1] {
            self.last_y_index.set(y_len - 1);
            y_len - 1
        } else {
            // Try accelerated search around last_y_index
            let start_y = self.last_y_index.get().saturating_sub(2);
            let mut found_y = false;
            let mut y_i = start_y;
            let max_y = y_len - 1;
            for k in 0..5 {
                let idx = (start_y + k).min(max_y);
                if idx == max_y {
                    break;
                }
                if self.y[idx] <= y_input && y_input <= self.y[idx + 1] {
                    y_i = idx;
                    found_y = true;
                    break;
                }
            }
            if found_y {
                y_i
            } else {
                self.find_bracket_binary(&self.y, y_input)
            }
        };

        self.last_x_index.set(x_idx);
        self.last_y_index.set(y_idx);

        let x1_idx = (x_idx + 1).min(x_len - 1);
        let y1_idx = (y_idx + 1).min(y_len - 1);

        let x0 = self.x[x_idx];
        let x1 = self.x[x1_idx];
        let y0 = self.y[y_idx];
        let y1 = self.y[y1_idx];

        let f00 = self.values[y_idx][x_idx];
        let f01 = self.values[y_idx][x1_idx];
        let f10 = self.values[y1_idx][x_idx];
        let f11 = self.values[y1_idx][x1_idx];

        if x0 == x1 || y0 == y1 {
            // Degenerate case, use nearest
            return f00;
        }

        // Bilinear interpolation
        let dx = ((x_input - x0).max(0.0).min(x1 - x0)) / (x1 - x0);
        let dy = ((y_input - y0).max(0.0).min(y1 - y0)) / (y1 - y0);

        (1.0 - dx) * (1.0 - dy) * f00
            + dx * (1.0 - dy) * f01
            + (1.0 - dx) * dy * f10
            + dx * dy * f11
    }

    /// Performs bilinear interpolation using binary search.
    /// Uses clamping for out-of-bounds values.
    pub fn interpolate_bilinear(&self, x_input: f64, y_input: f64) -> f64 {
        let x_len = self.x.len();
        let y_len = self.y.len();

        if x_len == 1 && y_len == 1 {
            return self.values[0][0];
        }

        // Find x indices
        let x0_idx = if x_input <= self.x[0] {
            0
        } else if x_input >= self.x[x_len - 1] {
            x_len - 1
        } else {
            self.find_bracket_binary(&self.x, x_input)
        };

        // Find y indices
        let y0_idx = if y_input <= self.y[0] {
            0
        } else if y_input >= self.y[y_len - 1] {
            y_len - 1
        } else {
            self.find_bracket_binary(&self.y, y_input)
        };

        let x1_idx = (x0_idx + 1).min(x_len - 1);
        let y1_idx = (y0_idx + 1).min(y_len - 1);

        let x0 = self.x[x0_idx];
        let x1 = self.x[x1_idx];
        let y0 = self.y[y0_idx];
        let y1 = self.y[y1_idx];

        let f00 = self.values[y0_idx][x0_idx];
        let f01 = self.values[y0_idx][x1_idx];
        let f10 = self.values[y1_idx][x0_idx];
        let f11 = self.values[y1_idx][x1_idx];

        if x0 == x1 || y0 == y1 {
            // Degenerate case, use nearest
            return f00;
        }

        // Bilinear interpolation with clamping
        let dx = if x1_idx == x0_idx {
            0.0
        } else {
            ((x_input - x0).max(0.0).min(x1 - x0)) / (x1 - x0)
        };
        let dy = if y1_idx == y0_idx {
            0.0
        } else {
            ((y_input - y0).max(0.0).min(y1 - y0)) / (y1 - y0)
        };

        (1.0 - dx) * (1.0 - dy) * f00
            + dx * (1.0 - dy) * f01
            + (1.0 - dx) * dy * f10
            + dx * dy * f11
    }

    /// Helper to find bracket index using binary search
    fn find_bracket_binary(&self, vec: &[f64], val: f64) -> usize {
        match vec.binary_search_by(|&v| v.partial_cmp(&val).unwrap_or(Ordering::Equal)) {
            Ok(idx) => idx,
            Err(insert_idx) => insert_idx.saturating_sub(1),
        }
    }
}

/// returns the sqare root of the input as is it was positive, but with the same sign as the input
pub fn signed_sqrt(x: f64) -> f64 {
    if x > 0. { x.sqrt() } else { -(-x).sqrt() }
}

#[cfg(test)]
mod tests {
    use super::*;
    use assert_float_eq::*;

    #[test]
    fn test_lookup_table_1d_single_point() {
        let table = LookupTable1D::new(vec![5.0], vec![10.0]);
        assert_float_relative_eq!(table.interpolate_linear_binary_search(5.0), 10.0);
        assert_float_relative_eq!(table.interpolate_linear_stateful(5.0), 10.0);
    }

    #[test]
    fn test_lookup_table_1d_clamping() {
        let table = LookupTable1D::new(vec![1.0, 2.0, 3.0], vec![10.0, 20.0, 30.0]);
        assert_float_relative_eq!(table.interpolate_linear_binary_search(0.5), 10.0); // Below min
        assert_float_relative_eq!(table.interpolate_linear_stateful(0.5), 10.0);
        assert_float_relative_eq!(table.interpolate_linear_linear_search(3.5), 30.0); // Above max
    }

    #[test]
    fn test_lookup_table_1d_interpolation() {
        let table = LookupTable1D::new(vec![1.0, 2.0, 3.0], vec![10.0, 20.0, 30.0]);
        assert_float_relative_eq!(table.interpolate_linear_linear_search(1.5), 15.0);
        assert_float_relative_eq!(table.interpolate_linear_binary_search(1.5), 15.0);
        assert_float_relative_eq!(table.interpolate_linear_stateful(1.5), 15.0);
    }

    #[test]
    fn test_lookup_table_2d_single_point() {
        let table = LookupTable2D::new(vec![1.0], vec![2.0], vec![vec![5.0]]);
        assert_float_relative_eq!(table.interpolate_bilinear(1.0, 2.0), 5.0);
    }

    #[test]
    fn test_lookup_table_2d_bilinear_interpolation() {
        let table = LookupTable2D::new(
            vec![1.0, 2.0, 3.0],
            vec![4.0, 5.0],
            vec![vec![10.0, 20.0, 30.0], vec![40.0, 50.0, 60.0]],
        );
        // Interpolation at (1.5, 4.5)
        // Expected: bilinear approx
        let result = table.interpolate_bilinear(1.5, 4.5);
        // Manual calc: (10+20)/2 at y=4, (40+50)/2 at y=5 -> (15,45), then interp: 15 + 0.5*(45-15) = 30
        assert_float_relative_eq!(result, 30.0, 1e-6);
    }

    #[test]
    fn test_lookup_table_2d_clamping() {
        let table = LookupTable2D::new(
            vec![1.0, 2.0],
            vec![1.0, 2.0],
            vec![vec![10.0, 20.0], vec![30.0, 40.0]],
        );
        assert_float_relative_eq!(table.interpolate_bilinear(0.5, 0.5), 10.0); // Below both
        assert_float_relative_eq!(table.interpolate_bilinear(2.5, 2.5), 40.0); // Above both
        assert_float_relative_eq!(table.interpolate_bilinear(0.5, 2.5), 30.0); // Mixed
    }

    #[test]
    fn test_lookup_table_2d_stateful_bilinear() {
        let table = LookupTable2D::new(
            vec![1.0, 2.0, 3.0],
            vec![4.0, 5.0],
            vec![vec![10.0, 20.0, 30.0], vec![40.0, 50.0, 60.0]],
        );
        // Test interpolation at (1.5, 4.5)
        let result_stateful = table.interpolate_bilinear_stateful(1.5, 4.5);
        let result_stateless = table.interpolate_bilinear(1.5, 4.5);
        assert_float_relative_eq!(result_stateful, result_stateless, 1e-6);
        assert_float_relative_eq!(result_stateful, 30.0, 1e-6); // Manual calc check

        // Test clamping
        assert_float_relative_eq!(table.interpolate_bilinear_stateful(0.5, 0.5), 10.0); // Below both
        assert_float_relative_eq!(table.interpolate_bilinear_stateful(3.5, 5.5), 60.0); // Above both
        assert_float_relative_eq!(table.interpolate_bilinear_stateful(0.5, 5.5), 40.0); // Mixed
    }

    #[test]
    fn test_signed_sqrt() {
        assert_float_relative_eq!(signed_sqrt(4.0), 2.0);
        assert_float_relative_eq!(signed_sqrt(-4.0), -2.0);
        assert_float_relative_eq!(signed_sqrt(-64.0), -8.0);
    }
}
