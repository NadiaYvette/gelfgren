//! Matrix operations

use crate::{Result, GelfgrenError};

/// A mathematical matrix stored in row-major order
#[derive(Debug, Clone, PartialEq)]
pub struct Matrix {
    data: Vec<f64>,
    rows: usize,
    cols: usize,
}

impl Matrix {
    /// Create a new matrix from a flat slice (row-major order)
    ///
    /// # Arguments
    /// * `data` - Flat array of matrix elements in row-major order
    /// * `rows` - Number of rows
    /// * `cols` - Number of columns
    pub fn new(data: &[f64], rows: usize, cols: usize) -> Result<Self> {
        if data.len() != rows * cols {
            return Err(GelfgrenError::InvalidDimension(
                format!("Data length {} does not match dimensions {}x{}", data.len(), rows, cols)
            ));
        }

        Ok(Self {
            data: data.to_vec(),
            rows,
            cols,
        })
    }

    /// Create a zero matrix
    pub fn zeros(rows: usize, cols: usize) -> Self {
        Self {
            data: vec![0.0; rows * cols],
            rows,
            cols,
        }
    }

    /// Create an identity matrix
    pub fn identity(size: usize) -> Self {
        let mut data = vec![0.0; size * size];
        for i in 0..size {
            data[i * size + i] = 1.0;
        }
        Self {
            data,
            rows: size,
            cols: size,
        }
    }

    /// Get the number of rows
    pub fn rows(&self) -> usize {
        self.rows
    }

    /// Get the number of columns
    pub fn cols(&self) -> usize {
        self.cols
    }

    /// Get an element at (row, col)
    pub fn get(&self, row: usize, col: usize) -> Result<f64> {
        if row >= self.rows || col >= self.cols {
            return Err(GelfgrenError::InvalidArgument(
                format!("Index ({}, {}) out of bounds for {}x{} matrix", row, col, self.rows, self.cols)
            ));
        }
        Ok(self.data[row * self.cols + col])
    }

    /// Set an element at (row, col)
    pub fn set(&mut self, row: usize, col: usize, value: f64) -> Result<()> {
        if row >= self.rows || col >= self.cols {
            return Err(GelfgrenError::InvalidArgument(
                format!("Index ({}, {}) out of bounds for {}x{} matrix", row, col, self.rows, self.cols)
            ));
        }
        self.data[row * self.cols + col] = value;
        Ok(())
    }

    /// Get a reference to the underlying data
    pub fn as_slice(&self) -> &[f64] {
        &self.data
    }

    /// Matrix multiplication
    pub fn multiply(&self, other: &Matrix) -> Result<Matrix> {
        if self.cols != other.rows {
            return Err(GelfgrenError::InvalidDimension(
                format!("Cannot multiply {}x{} by {}x{} matrix",
                    self.rows, self.cols, other.rows, other.cols)
            ));
        }

        let mut result = Matrix::zeros(self.rows, other.cols);

        for i in 0..self.rows {
            for j in 0..other.cols {
                let mut sum = 0.0;
                for k in 0..self.cols {
                    sum += self.get(i, k)? * other.get(k, j)?;
                }
                result.set(i, j, sum)?;
            }
        }

        Ok(result)
    }

    /// Transpose the matrix
    pub fn transpose(&self) -> Matrix {
        let mut result = Matrix::zeros(self.cols, self.rows);

        for i in 0..self.rows {
            for j in 0..self.cols {
                let val = self.get(i, j).unwrap();
                result.set(j, i, val).unwrap();
            }
        }

        result
    }

    /// Add two matrices element-wise
    pub fn add(&self, other: &Matrix) -> Result<Matrix> {
        if self.rows != other.rows || self.cols != other.cols {
            return Err(GelfgrenError::InvalidDimension(
                format!("Matrix dimensions must match: {}x{} != {}x{}",
                    self.rows, self.cols, other.rows, other.cols)
            ));
        }

        Ok(Matrix {
            data: self.data.iter()
                .zip(other.data.iter())
                .map(|(a, b)| a + b)
                .collect(),
            rows: self.rows,
            cols: self.cols,
        })
    }

    /// Scale matrix by a scalar
    pub fn scale(&self, scalar: f64) -> Matrix {
        Matrix {
            data: self.data.iter().map(|x| x * scalar).collect(),
            rows: self.rows,
            cols: self.cols,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_matrix_creation() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let m = Matrix::new(&data, 2, 3).unwrap();
        assert_eq!(m.rows(), 2);
        assert_eq!(m.cols(), 3);
        assert_eq!(m.get(0, 0).unwrap(), 1.0);
        assert_eq!(m.get(1, 2).unwrap(), 6.0);
    }

    #[test]
    fn test_identity_matrix() {
        let m = Matrix::identity(3);
        assert_eq!(m.get(0, 0).unwrap(), 1.0);
        assert_eq!(m.get(1, 1).unwrap(), 1.0);
        assert_eq!(m.get(2, 2).unwrap(), 1.0);
        assert_eq!(m.get(0, 1).unwrap(), 0.0);
    }

    #[test]
    fn test_matrix_multiply() {
        let a = Matrix::new(&[1.0, 2.0, 3.0, 4.0], 2, 2).unwrap();
        let b = Matrix::new(&[5.0, 6.0, 7.0, 8.0], 2, 2).unwrap();
        let c = a.multiply(&b).unwrap();

        // [1 2] [5 6]   [19 22]
        // [3 4] [7 8] = [43 50]
        assert_relative_eq!(c.get(0, 0).unwrap(), 19.0);
        assert_relative_eq!(c.get(0, 1).unwrap(), 22.0);
        assert_relative_eq!(c.get(1, 0).unwrap(), 43.0);
        assert_relative_eq!(c.get(1, 1).unwrap(), 50.0);
    }

    #[test]
    fn test_matrix_transpose() {
        let m = Matrix::new(&[1.0, 2.0, 3.0, 4.0, 5.0, 6.0], 2, 3).unwrap();
        let mt = m.transpose();
        assert_eq!(mt.rows(), 3);
        assert_eq!(mt.cols(), 2);
        assert_eq!(mt.get(0, 0).unwrap(), 1.0);
        assert_eq!(mt.get(2, 1).unwrap(), 6.0);
    }
}
