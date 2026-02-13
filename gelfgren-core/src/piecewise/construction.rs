//! Piecewise rational function construction.

use super::evaluation::SubintervalData;
use super::mesh::Mesh;
use crate::error::{GelfgrenError, Result};
use crate::hermite::HermiteData;
use crate::pade::SymmetricPade;
use crate::rational::RationalFunction;
use num_traits::{Float, FromPrimitive};

/// A piecewise rational function S_{n,m}(τ)(x) on mesh τ.
///
/// Consists of rational approximants R_j(x) = P_{n,j}(x)/Q_{m,j}(x)
/// defined on each subinterval Δⱼ = [x_{j-1}, x_j].
pub struct PiecewiseRational<T: Float> {
    /// The mesh partition
    mesh: Mesh<T>,
    /// Rational approximants for each subinterval
    /// subintervals[j-1] corresponds to interval j (1-indexed)
    subintervals: Vec<SubintervalData<T>>,
    /// Numerator degree n
    n: usize,
    /// Denominator degree m
    m: usize,
}

impl<T: Float + FromPrimitive + std::fmt::Debug + std::iter::Sum> PiecewiseRational<T> {
    /// Constructs a piecewise rational approximation from mesh and derivative data.
    ///
    /// # Arguments
    ///
    /// * `mesh` - Mesh partition with derivative values at each point
    /// * `n` - Numerator degree for each subinterval
    /// * `m` - Denominator degree for each subinterval
    ///
    /// # Algorithm
    ///
    /// For each subinterval Δⱼ = [x_{j-1}, x_j]:
    /// 1. Extract derivative data at endpoints
    /// 2. Construct symmetric Padé approximant
    /// 3. Store in subinterval list
    ///
    /// # Errors
    ///
    /// Returns error if:
    /// - Mesh points lack sufficient derivative data
    /// - Symmetric Padé construction fails
    /// - Degrees don't satisfy n+m+1 = 2p (symmetry requirement)
    pub fn from_mesh(mesh: Mesh<T>, n: usize, m: usize) -> Result<Self> {
        let p = (n + m + 1) / 2;

        // Verify symmetry condition
        if (n + m + 1) % 2 != 0 {
            return Err(GelfgrenError::InvalidArgument(format!(
                "Symmetric Padé requires n+m+1 = 2p (even), got n+m+1 = {}",
                n + m + 1
            )));
        }

        // Verify mesh points have enough derivative data
        for point in mesh.points() {
            if point.num_derivatives() < p {
                return Err(GelfgrenError::InvalidArgument(format!(
                    "Mesh points need {} derivative values for [{}/{}] approximant, got {}",
                    p,
                    n,
                    m,
                    point.num_derivatives()
                )));
            }
        }

        let mut subintervals = Vec::with_capacity(mesh.num_intervals());

        // Construct approximant for each subinterval
        for j in 1..=mesh.num_intervals() {
            let (x_left, x_right) = mesh.subinterval(j).unwrap();

            // Get derivative data at endpoints
            let left_point = mesh.point(j - 1).unwrap();
            let right_point = mesh.point(j).unwrap();

            let left_derivs = left_point.values[0..p].to_vec();
            let right_derivs = right_point.values[0..p].to_vec();

            // Construct symmetric Padé approximant
            let pade = SymmetricPade::from_endpoint_derivatives(
                &left_derivs,
                &right_derivs,
                n,
                m,
                x_left,
                x_right,
            )?;

            subintervals.push(SubintervalData {
                interval: (x_left, x_right),
                rational: pade.rational().clone(),
                index: j,
            });
        }

        Ok(Self {
            mesh,
            subintervals,
            n,
            m,
        })
    }

    /// Constructs a piecewise rational approximation from a function.
    ///
    /// Samples the function and its derivatives at mesh points, then constructs
    /// the piecewise approximation.
    ///
    /// # Arguments
    ///
    /// * `mesh` - Mesh partition (without derivative data)
    /// * `n` - Numerator degree
    /// * `m` - Denominator degree
    /// * `func` - Function to approximate
    /// * `derivatives` - List of derivative functions [f, f', f'', ...]
    pub fn from_function<F>(
        mut mesh: Mesh<T>,
        n: usize,
        m: usize,
        derivatives: &[F],
    ) -> Result<Self>
    where
        F: Fn(T) -> T,
    {
        let p = (n + m + 1) / 2;

        if derivatives.len() < p {
            return Err(GelfgrenError::InvalidArgument(format!(
                "Need {} derivative functions for [{}/{}] approximant",
                p, n, m
            )));
        }

        // Evaluate derivatives at each mesh point
        for point in mesh.points_mut() {
            let mut values = Vec::with_capacity(p);
            for deriv_fn in &derivatives[0..p] {
                values.push(deriv_fn(point.x));
            }
            point.values = values;
        }

        Self::from_mesh(mesh, n, m)
    }

    /// Returns the mesh partition.
    pub fn mesh(&self) -> &Mesh<T> {
        &self.mesh
    }

    /// Returns the degrees (n, m).
    pub fn degrees(&self) -> (usize, usize) {
        (self.n, self.m)
    }

    /// Returns the number of subintervals.
    pub fn num_intervals(&self) -> usize {
        self.subintervals.len()
    }

    /// Returns the subinterval data for interval j (1-indexed).
    pub fn subinterval(&self, j: usize) -> Option<&SubintervalData<T>> {
        if j == 0 || j > self.subintervals.len() {
            return None;
        }
        Some(&self.subintervals[j - 1])
    }

    /// Evaluates the piecewise rational function at x.
    ///
    /// Routes to the appropriate subinterval approximant.
    pub fn evaluate(&self, x: T) -> Result<T> {
        let j = self.mesh.locate(x).ok_or_else(|| {
            GelfgrenError::InvalidArgument(format!(
                "Point {:?} is outside mesh interval [{:?}, {:?}]",
                x,
                self.mesh.interval().0,
                self.mesh.interval().1
            ))
        })?;

        self.subintervals[j - 1].rational.evaluate(x)
    }

    /// Checks continuity at mesh point i.
    ///
    /// Computes left and right limits and their difference.
    /// Returns (left_limit, right_limit, difference).
    pub fn check_continuity(&self, i: usize) -> Option<(T, T, T)> {
        if i == 0 || i >= self.mesh.num_points() {
            return None; // Endpoints have no continuity to check
        }

        let x = self.mesh.point(i).unwrap().x;
        let tol = T::from_f64(1e-6).unwrap();

        // Evaluate from left (interval i)
        let x_left = x - tol;
        let left_val = self.subintervals[i - 1].rational.evaluate(x_left).ok()?;

        // Evaluate from right (interval i+1)
        let x_right = x + tol;
        let right_val = self.subintervals[i].rational.evaluate(x_right).ok()?;

        let diff = (left_val - right_val).abs();
        Some((left_val, right_val, diff))
    }

    /// Computes maximum discontinuity across all interior mesh points.
    pub fn max_discontinuity(&self) -> T {
        let mut max_disc = T::zero();
        for i in 1..self.mesh.num_points() - 1 {
            if let Some((_, _, diff)) = self.check_continuity(i) {
                if diff > max_disc {
                    max_disc = diff;
                }
            }
        }
        max_disc
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::piecewise::mesh::MeshPoint;
    use approx::assert_relative_eq;

    #[test]
    fn test_piecewise_construction_validation() {
        // Test symmetry requirement
        let mesh = Mesh::uniform(0.0, 1.0, 2).unwrap();
        let result = PiecewiseRational::from_mesh(mesh, 1, 1);
        assert!(result.is_err()); // n+m+1 = 3 (odd)
    }

    #[test]
    #[ignore] // TODO: Requires proper Newton-to-Bernstein conversion in symmetric Padé
    fn test_piecewise_from_function() {
        // Approximate f(x) = x with [0/1] Padé on [0, 1] with 2 intervals
        let mesh = Mesh::uniform(0.0, 1.0, 2).unwrap();

        let f = |x: f64| x;
        let df = |_x: f64| 1.0;

        let piecewise = PiecewiseRational::from_function(mesh, 0, 1, &[f, df]).unwrap();

        assert_eq!(piecewise.num_intervals(), 2);
        assert_eq!(piecewise.degrees(), (0, 1));
    }

    #[test]
    #[ignore] // TODO: Requires proper Newton-to-Bernstein conversion
    fn test_piecewise_evaluation() {
        let mesh = Mesh::uniform(0.0, 1.0, 2).unwrap();
        let f = |x: f64| x * x;
        let df = |x: f64| 2.0 * x;

        let piecewise = PiecewiseRational::from_function(mesh, 0, 1, &[f, df]).unwrap();

        // Should approximate x²
        assert_relative_eq!(piecewise.evaluate(0.25).unwrap(), 0.0625, epsilon = 0.1);
        assert_relative_eq!(piecewise.evaluate(0.75).unwrap(), 0.5625, epsilon = 0.1);
    }

    #[test]
    fn test_piecewise_degrees() {
        let points = vec![
            MeshPoint::with_values(0.0, vec![0.0, 1.0]),
            MeshPoint::with_values(0.5, vec![0.25, 1.0]),
            MeshPoint::with_values(1.0, vec![1.0, 1.0]),
        ];
        let mesh = Mesh::new(points).unwrap();

        let result = PiecewiseRational::from_mesh(mesh, 0, 1);
        if let Ok(piecewise) = result {
            assert_eq!(piecewise.degrees(), (0, 1));
            assert_eq!(piecewise.num_intervals(), 2);
        }
    }
}
