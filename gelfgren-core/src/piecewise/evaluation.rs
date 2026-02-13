//! Subinterval data and evaluation utilities.

use crate::rational::RationalFunction;
use num_traits::Float;

/// Data for a single subinterval of the piecewise function.
#[derive(Debug, Clone)]
pub struct SubintervalData<T: Float> {
    /// The interval [x_{j-1}, x_j]
    pub interval: (T, T),
    /// The rational approximant R_j(x) on this interval
    pub rational: RationalFunction<T>,
    /// Subinterval index (1-indexed, matching Gelfgren's notation)
    pub index: usize,
}

impl<T: Float> SubintervalData<T> {
    /// Returns the interval width h_j = x_j - x_{j-1}.
    pub fn width(&self) -> T {
        self.interval.1 - self.interval.0
    }

    /// Returns the midpoint (x_{j-1} + x_j)/2.
    pub fn midpoint(&self) -> T {
        (self.interval.0 + self.interval.1) / T::from(2.0).unwrap()
    }

    /// Checks if x is in this subinterval.
    pub fn contains(&self, x: T) -> bool {
        x >= self.interval.0 && x <= self.interval.1
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bernstein::BernsteinPolynomial;

    #[test]
    fn test_subinterval_data() {
        let p = BernsteinPolynomial::constant(1.0, 0.0, 1.0).unwrap();
        let q = BernsteinPolynomial::constant(1.0, 0.0, 1.0).unwrap();
        let rational = RationalFunction::new(p, q).unwrap();

        let data = SubintervalData {
            interval: (0.0, 1.0),
            rational,
            index: 1,
        };

        assert_eq!(data.width(), 1.0);
        assert_eq!(data.midpoint(), 0.5);
        assert!(data.contains(0.5));
        assert!(!data.contains(1.5));
    }
}
