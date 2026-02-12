//! Vector operations

use crate::{Result, GelfgrenError};

/// A mathematical vector
#[derive(Debug, Clone, PartialEq)]
pub struct Vector {
    data: Vec<f64>,
}

impl Vector {
    /// Create a new vector from a slice
    pub fn new(data: &[f64]) -> Self {
        Self {
            data: data.to_vec(),
        }
    }

    /// Create a zero vector of given length
    pub fn zeros(len: usize) -> Self {
        Self {
            data: vec![0.0; len],
        }
    }

    /// Create a vector of ones
    pub fn ones(len: usize) -> Self {
        Self {
            data: vec![1.0; len],
        }
    }

    /// Get the length of the vector
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// Check if the vector is empty
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Get a reference to the underlying data
    pub fn as_slice(&self) -> &[f64] {
        &self.data
    }

    /// Get a mutable reference to the underlying data
    pub fn as_mut_slice(&mut self) -> &mut [f64] {
        &mut self.data
    }

    /// Compute the dot product with another vector
    pub fn dot(&self, other: &Vector) -> Result<f64> {
        super::dot(&self.data, &other.data)
    }

    /// Compute the Euclidean norm
    pub fn norm(&self) -> f64 {
        super::norm(&self.data)
    }

    /// Add two vectors element-wise
    pub fn add(&self, other: &Vector) -> Result<Vector> {
        if self.len() != other.len() {
            return Err(GelfgrenError::InvalidDimension(
                format!("Vector lengths must match: {} != {}", self.len(), other.len())
            ));
        }

        Ok(Vector {
            data: self.data.iter()
                .zip(other.data.iter())
                .map(|(a, b)| a + b)
                .collect(),
        })
    }

    /// Subtract two vectors element-wise
    pub fn sub(&self, other: &Vector) -> Result<Vector> {
        if self.len() != other.len() {
            return Err(GelfgrenError::InvalidDimension(
                format!("Vector lengths must match: {} != {}", self.len(), other.len())
            ));
        }

        Ok(Vector {
            data: self.data.iter()
                .zip(other.data.iter())
                .map(|(a, b)| a - b)
                .collect(),
        })
    }

    /// Multiply vector by a scalar
    pub fn scale(&self, scalar: f64) -> Vector {
        Vector {
            data: self.data.iter().map(|x| x * scalar).collect(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_vector_creation() {
        let v = Vector::new(&[1.0, 2.0, 3.0]);
        assert_eq!(v.len(), 3);
        assert_eq!(v.as_slice(), &[1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_vector_zeros() {
        let v = Vector::zeros(5);
        assert_eq!(v.len(), 5);
        assert!(v.as_slice().iter().all(|&x| x == 0.0));
    }

    #[test]
    fn test_vector_add() {
        let v1 = Vector::new(&[1.0, 2.0, 3.0]);
        let v2 = Vector::new(&[4.0, 5.0, 6.0]);
        let result = v1.add(&v2).unwrap();
        assert_eq!(result.as_slice(), &[5.0, 7.0, 9.0]);
    }

    #[test]
    fn test_vector_scale() {
        let v = Vector::new(&[1.0, 2.0, 3.0]);
        let result = v.scale(2.0);
        assert_eq!(result.as_slice(), &[2.0, 4.0, 6.0]);
    }

    #[test]
    fn test_vector_norm() {
        let v = Vector::new(&[3.0, 4.0]);
        assert_relative_eq!(v.norm(), 5.0);
    }
}
