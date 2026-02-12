//! Statistical functions
//!
//! This module provides descriptive statistics, probability distributions,
//! and hypothesis testing functions.

use crate::Result;

/// Compute the mean (average) of a dataset
///
/// # Arguments
/// * `data` - Input data slice
///
/// # Returns
/// The arithmetic mean
///
/// # Errors
/// Returns `InvalidArgument` if the data slice is empty
pub fn mean(data: &[f64]) -> Result<f64> {
    if data.is_empty() {
        return Err(crate::GelfgrenError::InvalidArgument(
            "Cannot compute mean of empty dataset".to_string()
        ));
    }

    Ok(data.iter().sum::<f64>() / data.len() as f64)
}

/// Compute the variance of a dataset
///
/// # Arguments
/// * `data` - Input data slice
/// * `ddof` - Delta degrees of freedom (0 for population, 1 for sample)
///
/// # Returns
/// The variance
pub fn variance(data: &[f64], ddof: usize) -> Result<f64> {
    if data.is_empty() {
        return Err(crate::GelfgrenError::InvalidArgument(
            "Cannot compute variance of empty dataset".to_string()
        ));
    }

    if data.len() <= ddof {
        return Err(crate::GelfgrenError::InvalidArgument(
            format!("Dataset size {} must be greater than ddof {}", data.len(), ddof)
        ));
    }

    let m = mean(data)?;
    let sum_sq_diff: f64 = data.iter().map(|x| (x - m).powi(2)).sum();
    Ok(sum_sq_diff / (data.len() - ddof) as f64)
}

/// Compute the standard deviation of a dataset
///
/// # Arguments
/// * `data` - Input data slice
/// * `ddof` - Delta degrees of freedom (0 for population, 1 for sample)
///
/// # Returns
/// The standard deviation
pub fn std_dev(data: &[f64], ddof: usize) -> Result<f64> {
    Ok(variance(data, ddof)?.sqrt())
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_mean() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let result = mean(&data).unwrap();
        assert_relative_eq!(result, 3.0);
    }

    #[test]
    fn test_variance() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let result = variance(&data, 1).unwrap();  // Sample variance
        assert_relative_eq!(result, 2.5);
    }

    #[test]
    fn test_std_dev() {
        let data = vec![2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
        let result = std_dev(&data, 1).unwrap();
        assert_relative_eq!(result, 2.138, epsilon = 0.001);
    }
}
