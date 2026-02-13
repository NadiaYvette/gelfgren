//! Mesh partition management.

use crate::error::{GelfgrenError, Result};
use num_traits::{Float, FromPrimitive};

/// A mesh point with optional derivative information.
#[derive(Debug, Clone)]
pub struct MeshPoint<T: Float> {
    /// The x-coordinate of the mesh point
    pub x: T,
    /// Function and derivative values at this point: [f(x), f'(x), f''(x), ...]
    /// Optional - may be empty if only geometry is needed
    pub values: Vec<T>,
}

impl<T: Float> MeshPoint<T> {
    /// Creates a mesh point with position only.
    pub fn new(x: T) -> Self {
        Self {
            x,
            values: Vec::new(),
        }
    }

    /// Creates a mesh point with position and derivative values.
    pub fn with_values(x: T, values: Vec<T>) -> Self {
        Self { x, values }
    }

    /// Returns the function value f(x) if available.
    pub fn function_value(&self) -> Option<T> {
        self.values.first().copied()
    }

    /// Returns the derivative value f'(x) if available.
    pub fn derivative(&self, order: usize) -> Option<T> {
        self.values.get(order).copied()
    }

    /// Returns the number of available derivative values.
    pub fn num_derivatives(&self) -> usize {
        self.values.len()
    }
}

/// A mesh partition of an interval [a, b].
///
/// Represents τ: a = x₀ < x₁ < x₂ < ... < x_N = b
#[derive(Debug, Clone)]
pub struct Mesh<T: Float> {
    /// Mesh points in ascending order
    points: Vec<MeshPoint<T>>,
}

impl<T: Float + FromPrimitive> Mesh<T> {
    /// Creates a mesh from a list of mesh points.
    ///
    /// # Arguments
    ///
    /// * `points` - Mesh points (must be in ascending order)
    ///
    /// # Errors
    ///
    /// Returns error if:
    /// - Less than 2 points provided
    /// - Points are not strictly increasing
    pub fn new(points: Vec<MeshPoint<T>>) -> Result<Self> {
        if points.len() < 2 {
            return Err(GelfgrenError::InvalidArgument(
                "Mesh must have at least 2 points".to_string(),
            ));
        }

        // Verify strictly increasing
        for i in 1..points.len() {
            if points[i].x <= points[i - 1].x {
                return Err(GelfgrenError::InvalidArgument(
                    "Mesh points must be strictly increasing".to_string(),
                ));
            }
        }

        Ok(Self { points })
    }

    /// Creates a uniform mesh on [a, b] with n subintervals.
    ///
    /// Generates n+1 equally-spaced points: x_i = a + i·h where h = (b-a)/n
    pub fn uniform(a: T, b: T, n: usize) -> Result<Self>
    where
        T: std::fmt::Debug,
    {
        if n == 0 {
            return Err(GelfgrenError::InvalidArgument(
                "Must have at least 1 subinterval".to_string(),
            ));
        }

        if a >= b {
            return Err(GelfgrenError::InvalidArgument(format!(
                "Invalid interval: a={:?} must be less than b={:?}",
                a, b
            )));
        }

        let h = (b - a) / T::from_usize(n).unwrap();
        let mut points = Vec::with_capacity(n + 1);

        for i in 0..=n {
            let x = a + T::from_usize(i).unwrap() * h;
            points.push(MeshPoint::new(x));
        }

        Ok(Self { points })
    }

    /// Creates a Chebyshev mesh on [a, b] with n+1 points.
    ///
    /// Uses Chebyshev points of the second kind: x_i = cos(iπ/n)
    /// mapped to [a, b]. These points cluster near endpoints for
    /// better approximation of functions with boundary layers.
    pub fn chebyshev(a: T, b: T, n: usize) -> Result<Self>
    where
        T: std::fmt::Debug,
    {
        if n == 0 {
            return Err(GelfgrenError::InvalidArgument(
                "Must have at least 1 subinterval".to_string(),
            ));
        }

        let mid = (a + b) / T::from_usize(2).unwrap();
        let half_width = (b - a) / T::from_usize(2).unwrap();

        let mut points = Vec::with_capacity(n + 1);
        let pi = T::from_f64(std::f64::consts::PI).unwrap();
        let n_f = T::from_usize(n).unwrap();

        for i in 0..=n {
            let i_f = T::from_usize(i).unwrap();
            let theta = i_f * pi / n_f;
            let xi = -theta.cos(); // Map to [-1, 1]
            let x = mid + half_width * xi;
            points.push(MeshPoint::new(x));
        }

        // Sort to ensure ascending order (Chebyshev points go from b to a)
        points.sort_by(|p1, p2| p1.x.partial_cmp(&p2.x).unwrap());

        Ok(Self { points })
    }

    /// Returns the number of mesh points (N+1).
    #[inline]
    pub fn num_points(&self) -> usize {
        self.points.len()
    }

    /// Returns the number of subintervals (N).
    #[inline]
    pub fn num_intervals(&self) -> usize {
        self.points.len() - 1
    }

    /// Returns the mesh points.
    #[inline]
    pub fn points(&self) -> &[MeshPoint<T>] {
        &self.points
    }

    /// Returns mutable access to the mesh points.
    #[inline]
    pub fn points_mut(&mut self) -> &mut [MeshPoint<T>] {
        &mut self.points
    }

    /// Returns the i-th mesh point.
    #[inline]
    pub fn point(&self, i: usize) -> Option<&MeshPoint<T>> {
        self.points.get(i)
    }

    /// Returns the interval [a, b].
    pub fn interval(&self) -> (T, T) {
        (self.points.first().unwrap().x, self.points.last().unwrap().x)
    }

    /// Returns the j-th subinterval [x_{j-1}, x_j].
    pub fn subinterval(&self, j: usize) -> Option<(T, T)> {
        if j == 0 || j > self.num_intervals() {
            return None;
        }
        Some((self.points[j - 1].x, self.points[j].x))
    }

    /// Finds which subinterval contains x.
    ///
    /// Returns j such that x ∈ [x_{j-1}, x_j], or None if x is outside [a, b].
    pub fn locate(&self, x: T) -> Option<usize> {
        let (a, b) = self.interval();
        if x < a || x > b {
            return None;
        }

        // Linear search for the containing interval
        for j in 1..=self.num_intervals() {
            let (x_left, x_right) = self.subinterval(j).unwrap();
            if x >= x_left && x <= x_right {
                return Some(j);
            }
        }

        None
    }

    /// Computes the maximum mesh width h = max_j |x_j - x_{j-1}|.
    pub fn max_width(&self) -> T {
        let mut max_h = T::zero();
        for j in 1..self.num_points() {
            let h = self.points[j].x - self.points[j - 1].x;
            if h > max_h {
                max_h = h;
            }
        }
        max_h
    }

    /// Computes the minimum mesh width.
    pub fn min_width(&self) -> T {
        let mut min_h = T::infinity();
        for j in 1..self.num_points() {
            let h = self.points[j].x - self.points[j - 1].x;
            if h < min_h {
                min_h = h;
            }
        }
        min_h
    }

    /// Refines the mesh by subdividing each interval in half.
    ///
    /// Creates a new mesh with 2N+1 points.
    pub fn refine_uniform(&self) -> Self {
        let mut new_points = Vec::with_capacity(2 * self.num_points() - 1);

        for j in 0..self.num_intervals() {
            let (x_left, x_right) = self.subinterval(j + 1).unwrap();
            let x_mid = (x_left + x_right) / T::from_usize(2).unwrap();

            new_points.push(self.points[j].clone());
            new_points.push(MeshPoint::new(x_mid));
        }
        new_points.push(self.points.last().unwrap().clone());

        Self { points: new_points }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_mesh_point_creation() {
        let p = MeshPoint::new(1.0);
        assert_eq!(p.x, 1.0);
        assert_eq!(p.num_derivatives(), 0);

        let p = MeshPoint::with_values(2.0, vec![3.0, 4.0]);
        assert_eq!(p.x, 2.0);
        assert_eq!(p.function_value(), Some(3.0));
        assert_eq!(p.derivative(0), Some(3.0));
        assert_eq!(p.derivative(1), Some(4.0));
    }

    #[test]
    fn test_uniform_mesh() {
        let mesh = Mesh::uniform(0.0, 1.0, 4).unwrap();
        assert_eq!(mesh.num_points(), 5);
        assert_eq!(mesh.num_intervals(), 4);

        for i in 0..5 {
            assert_relative_eq!(mesh.point(i).unwrap().x, 0.25 * i as f64, epsilon = 1e-10);
        }

        assert_relative_eq!(mesh.max_width(), 0.25, epsilon = 1e-10);
    }

    #[test]
    fn test_chebyshev_mesh() {
        let mesh = Mesh::chebyshev(0.0, 1.0, 4).unwrap();
        assert_eq!(mesh.num_points(), 5);

        // Verify points are in ascending order
        for i in 1..5 {
            assert!(mesh.point(i).unwrap().x > mesh.point(i - 1).unwrap().x);
        }

        // Endpoints should be close to 0 and 1
        assert_relative_eq!(mesh.point(0).unwrap().x, 0.0, epsilon = 1e-10);
        assert_relative_eq!(mesh.point(4).unwrap().x, 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_subinterval_access() {
        let mesh = Mesh::uniform(0.0, 1.0, 2).unwrap();

        assert_eq!(mesh.subinterval(1), Some((0.0, 0.5)));
        assert_eq!(mesh.subinterval(2), Some((0.5, 1.0)));
        assert_eq!(mesh.subinterval(0), None);
        assert_eq!(mesh.subinterval(3), None);
    }

    #[test]
    fn test_locate() {
        let mesh = Mesh::uniform(0.0, 1.0, 4).unwrap();

        assert_eq!(mesh.locate(0.1), Some(1));
        assert_eq!(mesh.locate(0.3), Some(2));
        assert_eq!(mesh.locate(0.5), Some(2)); // On boundary
        assert_eq!(mesh.locate(0.9), Some(4));
        assert_eq!(mesh.locate(-0.1), None);
        assert_eq!(mesh.locate(1.1), None);
    }

    #[test]
    fn test_refine_uniform() {
        let mesh = Mesh::uniform(0.0, 1.0, 2).unwrap();
        let refined = mesh.refine_uniform();

        assert_eq!(refined.num_intervals(), 4);
        assert_relative_eq!(refined.max_width(), 0.25, epsilon = 1e-10);
    }

    #[test]
    fn test_mesh_validation() {
        // Non-increasing points
        let points = vec![MeshPoint::new(1.0), MeshPoint::new(0.0)];
        assert!(Mesh::new(points).is_err());

        // Too few points
        let points = vec![MeshPoint::new(0.0)];
        assert!(Mesh::new(points).is_err());
    }

    #[test]
    fn test_min_max_width() {
        let points = vec![
            MeshPoint::new(0.0),
            MeshPoint::new(0.1),
            MeshPoint::new(0.5),
            MeshPoint::new(1.0),
        ];
        let mesh = Mesh::new(points).unwrap();

        assert_relative_eq!(mesh.min_width(), 0.1, epsilon = 1e-10);
        assert_relative_eq!(mesh.max_width(), 0.5, epsilon = 1e-10);
    }
}
