//! Python bindings for Gelfgren library using PyO3.
//!
//! Provides Pythonic interface with NumPy integration and exception handling.

use gelfgren_core::bernstein::BernsteinPolynomial as CoreBernstein;
use gelfgren_core::rational::RationalFunction as CoreRational;
use gelfgren_core::pade::{PadeApproximant as CorePade, TwoPointPade as CoreTwoPointPade, HermiteConstraints as CoreHermiteConstraints};
use gelfgren_core::piecewise::{Mesh as CoreMesh, MeshPoint, PiecewiseRational as CorePiecewiseRational};
use gelfgren_core::error::GelfgrenError;

use pyo3::prelude::*;
use pyo3::exceptions::{PyValueError, PyRuntimeError};
use numpy::{PyArray1, PyArray2, PyReadonlyArray1, ToPyArray};

/// Convert GelfgrenError to Python exception
fn to_py_err(err: GelfgrenError) -> PyErr {
    match err {
        GelfgrenError::InvalidArgument(msg) => PyValueError::new_err(msg),
        GelfgrenError::InvalidDimension(msg) => PyValueError::new_err(msg),
        GelfgrenError::DivisionByZero => PyRuntimeError::new_err("Division by zero"),
        GelfgrenError::SingularMatrix => PyRuntimeError::new_err("Singular matrix"),
        GelfgrenError::NotImplemented(msg) => PyRuntimeError::new_err(format!("Not implemented: {}", msg)),
        GelfgrenError::ConvergenceFailure(msg) => PyRuntimeError::new_err(format!("Convergence failure: {}", msg)),
        GelfgrenError::Overflow => PyRuntimeError::new_err("Numerical overflow"),
        GelfgrenError::Underflow => PyRuntimeError::new_err("Numerical underflow"),
    }
}

/// Bernstein polynomial in Python.
///
/// Provides numerically stable polynomial representation using Bernstein basis.
///
/// Parameters
/// ----------
/// coeffs : array_like
///     Unscaled Bernstein coefficients [C₀, C₁, ..., Cₙ]
/// a : float
///     Left endpoint of interval
/// b : float
///     Right endpoint of interval
///
/// Examples
/// --------
/// >>> import gelfgren
/// >>> import numpy as np
/// >>> # Create polynomial P(x) = 1 + 2x + 3x²
/// >>> p = gelfgren.BernsteinPolynomial([1.0, 2.0, 3.0], 0.0, 1.0)
/// >>> p.evaluate(0.5)
/// 1.5
/// >>> p.degree()
/// 2
#[pyclass(name = "BernsteinPolynomial")]
struct PyBernstein {
    inner: CoreBernstein<f64>,
}

#[pymethods]
impl PyBernstein {
    #[new]
    fn new(coeffs: PyReadonlyArray1<f64>, a: f64, b: f64) -> PyResult<Self> {
        let coeffs_vec = coeffs.as_slice()?.to_vec();
        let inner = CoreBernstein::new(coeffs_vec, a, b)
            .map_err(to_py_err)?;
        Ok(Self { inner })
    }

    /// Evaluate the polynomial at a point or array of points.
    ///
    /// Parameters
    /// ----------
    /// x : float or array_like
    ///     Evaluation point(s)
    ///
    /// Returns
    /// -------
    /// float or ndarray
    ///     P(x)
    fn evaluate<'py>(&self, py: Python<'py>, x: PyReadonlyArray1<f64>) -> Bound<'py, PyArray1<f64>> {
        let x_slice = x.as_slice().unwrap();
        let result: Vec<f64> = x_slice.iter()
            .map(|&xi| self.inner.evaluate(xi))
            .collect();
        result.to_pyarray_bound(py)
    }

    /// Evaluate at a single point (convenience method).
    fn eval_scalar(&self, x: f64) -> f64 {
        self.inner.evaluate(x)
    }

    /// Get the degree of the polynomial.
    fn degree(&self) -> usize {
        self.inner.degree()
    }

    /// Get the interval [a, b].
    fn interval(&self) -> (f64, f64) {
        self.inner.interval()
    }

    /// Compute the derivative polynomial.
    fn derivative(&self) -> PyResult<Self> {
        let inner = self.inner.derivative().map_err(to_py_err)?;
        Ok(Self { inner })
    }

    /// Compute the integral (antiderivative) polynomial.
    fn integral(&self) -> PyResult<Self> {
        let inner = self.inner.integral().map_err(to_py_err)?;
        Ok(Self { inner })
    }

    /// Elevate the polynomial degree by 1.
    fn elevate(&self) -> Self {
        let inner = self.inner.degree_elevate();
        Self { inner }
    }

    /// Add two polynomials.
    fn __add__(&self, other: &Self) -> PyResult<Self> {
        let inner = (self.inner.clone() + other.inner.clone())
            .map_err(to_py_err)?;
        Ok(Self { inner })
    }

    /// Subtract two polynomials.
    fn __sub__(&self, other: &Self) -> PyResult<Self> {
        let inner = (self.inner.clone() - other.inner.clone())
            .map_err(to_py_err)?;
        Ok(Self { inner })
    }

    /// Multiply two polynomials.
    fn __mul__(&self, other: &Self) -> PyResult<Self> {
        let inner = (self.inner.clone() * other.inner.clone())
            .map_err(to_py_err)?;
        Ok(Self { inner })
    }

    /// String representation.
    fn __repr__(&self) -> String {
        let (a, b) = self.inner.interval();
        format!("BernsteinPolynomial(degree={}, interval=[{}, {}])",
                self.inner.degree(), a, b)
    }
}

/// Rational function P(x)/Q(x) in Python.
///
/// Parameters
/// ----------
/// numerator : BernsteinPolynomial
///     Numerator polynomial P(x)
/// denominator : BernsteinPolynomial
///     Denominator polynomial Q(x)
///
/// Examples
/// --------
/// >>> num = gelfgren.BernsteinPolynomial([1.0, 1.0], 0.0, 1.0)
/// >>> den = gelfgren.BernsteinPolynomial([1.0, 2.0], 0.0, 1.0)
/// >>> r = gelfgren.RationalFunction(num, den)
/// >>> r.evaluate(0.5)
/// 0.75
#[pyclass(name = "RationalFunction")]
struct PyRational {
    inner: CoreRational<f64>,
}

#[pymethods]
impl PyRational {
    #[new]
    fn new(numerator: &PyBernstein, denominator: &PyBernstein) -> PyResult<Self> {
        let inner = CoreRational::new(
            numerator.inner.clone(),
            denominator.inner.clone(),
        ).map_err(to_py_err)?;
        Ok(Self { inner })
    }

    /// Evaluate the rational function at a point.
    ///
    /// Raises
    /// ------
    /// RuntimeError
    ///     If denominator is zero at the evaluation point
    fn evaluate<'py>(&self, py: Python<'py>, x: PyReadonlyArray1<f64>) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let x_slice = x.as_slice().unwrap();
        let result: Result<Vec<f64>, _> = x_slice.iter()
            .map(|&xi| self.inner.evaluate(xi))
            .collect();
        let values = result.map_err(to_py_err)?;
        Ok(values.to_pyarray_bound(py))
    }

    /// Evaluate at a single point (convenience method).
    fn eval_scalar(&self, x: f64) -> PyResult<f64> {
        self.inner.evaluate(x).map_err(to_py_err)
    }

    /// Compute the derivative using the quotient rule.
    fn derivative(&self) -> PyResult<Self> {
        let inner = self.inner.derivative().map_err(to_py_err)?;
        Ok(Self { inner })
    }

    fn __repr__(&self) -> String {
        "RationalFunction(P/Q)".to_string()
    }
}

/// Padé approximant [n/m] in Python.
///
/// Parameters
/// ----------
/// coeffs : array_like
///     Power series coefficients [c₀, c₁, ..., c_{n+m}]
/// n : int
///     Numerator degree
/// m : int
///     Denominator degree
/// center : float
///     Expansion center
/// a : float
///     Left endpoint of interval
/// b : float
///     Right endpoint of interval
///
/// Examples
/// --------
/// >>> # Approximate exp(x) with [2/2] Padé
/// >>> coeffs = [1.0, 1.0, 0.5, 1.0/6.0, 1.0/24.0]
/// >>> pade = gelfgren.PadeApproximant(coeffs, 2, 2, 0.0, -1.0, 1.0)
/// >>> pade.evaluate(0.5)
/// 1.648...
#[pyclass(name = "PadeApproximant")]
struct PyPade {
    inner: CorePade<f64>,
}

#[pymethods]
impl PyPade {
    #[new]
    fn new(coeffs: PyReadonlyArray1<f64>, n: usize, m: usize,
           center: f64, a: f64, b: f64) -> PyResult<Self> {
        let coeffs_vec = coeffs.as_slice()?.to_vec();
        let inner = CorePade::from_series(&coeffs_vec, n, m, center, (a, b))
            .map_err(to_py_err)?;
        Ok(Self { inner })
    }

    /// Evaluate the Padé approximant at a point.
    fn evaluate<'py>(&self, py: Python<'py>, x: PyReadonlyArray1<f64>) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let x_slice = x.as_slice().unwrap();
        let result: Result<Vec<f64>, _> = x_slice.iter()
            .map(|&xi| self.inner.evaluate(xi))
            .collect();
        let values = result.map_err(to_py_err)?;
        Ok(values.to_pyarray_bound(py))
    }

    fn eval_scalar(&self, x: f64) -> PyResult<f64> {
        self.inner.evaluate(x).map_err(to_py_err)
    }

    fn __repr__(&self) -> String {
        "PadeApproximant([n/m])".to_string()
    }
}

/// Mesh partition for piecewise approximation.
///
/// Examples
/// --------
/// >>> mesh = gelfgren.Mesh.uniform(0.0, 1.0, 4)
/// >>> mesh = gelfgren.Mesh.chebyshev(0.0, 1.0, 4)
#[pyclass(name = "Mesh")]
#[derive(Clone)]
struct PyMesh {
    inner: CoreMesh<f64>,
}

#[pymethods]
impl PyMesh {
    /// Create a uniform mesh with equal spacing.
    #[staticmethod]
    fn uniform(a: f64, b: f64, n: usize) -> PyResult<Self> {
        let inner = CoreMesh::uniform(a, b, n).map_err(to_py_err)?;
        Ok(Self { inner })
    }

    /// Create a Chebyshev mesh with boundary clustering.
    #[staticmethod]
    fn chebyshev(a: f64, b: f64, n: usize) -> PyResult<Self> {
        let inner = CoreMesh::chebyshev(a, b, n).map_err(to_py_err)?;
        Ok(Self { inner })
    }

    /// Sample a function and its derivatives at mesh points.
    ///
    /// Parameters
    /// ----------
    /// funcs : list of callable
    ///     List of functions [f, f', f'', ...] to evaluate at each mesh point
    ///
    /// Returns
    /// -------
    /// Mesh
    ///     New mesh with derivative values populated
    fn sample_functions(&self, funcs: Vec<PyObject>) -> PyResult<Self> {
        Python::with_gil(|py| {
            let mut new_mesh = self.inner.clone();
            let num_derivs = funcs.len();

            for point in new_mesh.points_mut() {
                let mut values = Vec::with_capacity(num_derivs);
                for func in &funcs {
                    let result: f64 = func.call1(py, (point.x,))?.extract(py)?;
                    values.push(result);
                }
                point.values = values;
            }

            Ok(Self { inner: new_mesh })
        })
    }

    fn __repr__(&self) -> String {
        let (a, b) = self.inner.interval();
        format!("Mesh(interval=[{}, {}], points={})", a, b, self.inner.num_points())
    }
}

/// Two-point Padé approximant.
///
/// Constructs rational approximant matching derivatives at both endpoints.
///
/// Examples
/// --------
/// >>> import numpy as np
/// >>> # Approximate f(x) = exp(x) on [0, 1] with [2/1] rational
/// >>> left_derivs = [1.0, 1.0]  # f(0), f'(0)
/// >>> right_derivs = [np.e, np.e]  # f(1), f'(1)
/// >>> pade = gelfgren.TwoPointPade.from_derivatives(
/// ...     left_derivs, right_derivs, 2, 1, 0.0, 1.0)
/// >>> pade.evaluate(0.5)
#[pyclass(name = "TwoPointPade")]
struct PyTwoPointPade {
    inner: CoreTwoPointPade<f64>,
}

#[pymethods]
impl PyTwoPointPade {
    /// Construct from endpoint derivatives.
    ///
    /// Parameters
    /// ----------
    /// left_derivs : array_like
    ///     Derivatives at left endpoint [f(a), f'(a), ...]
    /// right_derivs : array_like
    ///     Derivatives at right endpoint [f(b), f'(b), ...]
    /// n : int
    ///     Numerator degree
    /// m : int
    ///     Denominator degree
    /// a : float
    ///     Left endpoint
    /// b : float
    ///     Right endpoint
    ///
    /// Notes
    /// -----
    /// Must satisfy n + m + 1 = 2p (even) where p = len(left_derivs)
    #[staticmethod]
    fn from_derivatives(
        left_derivs: PyReadonlyArray1<f64>,
        right_derivs: PyReadonlyArray1<f64>,
        n: usize,
        m: usize,
        a: f64,
        b: f64,
    ) -> PyResult<Self> {
        let left_vec = left_derivs.as_slice()?.to_vec();
        let right_vec = right_derivs.as_slice()?.to_vec();

        let inner = CoreTwoPointPade::from_endpoint_derivatives(
            &left_vec,
            &right_vec,
            n,
            m,
            a,
            b,
        ).map_err(to_py_err)?;

        Ok(Self { inner })
    }

    /// Evaluate at a point or array of points.
    fn evaluate<'py>(&self, py: Python<'py>, x: PyReadonlyArray1<f64>) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let x_slice = x.as_slice().unwrap();
        let result: Result<Vec<f64>, _> = x_slice.iter()
            .map(|&xi| self.inner.evaluate(xi))
            .collect();
        let values = result.map_err(to_py_err)?;
        Ok(values.to_pyarray_bound(py))
    }

    /// Evaluate at a single point.
    fn eval_scalar(&self, x: f64) -> PyResult<f64> {
        self.inner.evaluate(x).map_err(to_py_err)
    }

    /// Get the underlying rational function.
    fn rational(&self) -> PyRational {
        PyRational {
            inner: self.inner.rational().clone(),
        }
    }

    fn __repr__(&self) -> String {
        "TwoPointPade([n/m])".to_string()
    }
}

/// Piecewise rational approximation.
///
/// A piecewise function where each subinterval uses a rational approximant.
///
/// Examples
/// --------
/// >>> import numpy as np
/// >>> mesh = gelfgren.Mesh.uniform(0.0, 1.0, 4)
/// >>> f = lambda x: np.exp(x)
/// >>> df = lambda x: np.exp(x)
/// >>> mesh = mesh.sample_functions([f, df])
/// >>> pw = gelfgren.PiecewiseRational.from_mesh(mesh, 2, 1)
/// >>> pw.evaluate(np.array([0.5]))
#[pyclass(name = "PiecewiseRational")]
struct PyPiecewiseRational {
    inner: CorePiecewiseRational<f64>,
}

#[pymethods]
impl PyPiecewiseRational {
    /// Construct from mesh with derivative data.
    ///
    /// Parameters
    /// ----------
    /// mesh : Mesh
    ///     Mesh partition with derivative values at each point
    /// n : int
    ///     Numerator degree for each subinterval
    /// m : int
    ///     Denominator degree for each subinterval
    ///
    /// Notes
    /// -----
    /// Must satisfy n + m + 1 = 2p (even) where p is the number of
    /// derivative values at each mesh point
    #[staticmethod]
    fn from_mesh(mesh: &PyMesh, n: usize, m: usize) -> PyResult<Self> {
        let inner = CorePiecewiseRational::from_mesh(mesh.inner.clone(), n, m)
            .map_err(to_py_err)?;
        Ok(Self { inner })
    }

    /// Construct from a function and its derivatives.
    ///
    /// Parameters
    /// ----------
    /// mesh : Mesh
    ///     Mesh partition (without derivative data)
    /// n : int
    ///     Numerator degree
    /// m : int
    ///     Denominator degree
    /// funcs : list of callable
    ///     List of functions [f, f', f'', ...] to sample at mesh points
    #[staticmethod]
    fn from_function(
        mesh: &PyMesh,
        n: usize,
        m: usize,
        funcs: Vec<PyObject>,
    ) -> PyResult<Self> {
        Python::with_gil(|py| {
            // Sample functions at mesh points
            let sampled_mesh = mesh.sample_functions(funcs)?;
            Self::from_mesh(&sampled_mesh, n, m)
        })
    }

    /// Evaluate at a point or array of points.
    fn evaluate<'py>(&self, py: Python<'py>, x: PyReadonlyArray1<f64>) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let x_slice = x.as_slice().unwrap();
        let result: Result<Vec<f64>, _> = x_slice.iter()
            .map(|&xi| self.inner.evaluate(xi))
            .collect();
        let values = result.map_err(to_py_err)?;
        Ok(values.to_pyarray_bound(py))
    }

    /// Evaluate at a single point.
    fn eval_scalar(&self, x: f64) -> PyResult<f64> {
        self.inner.evaluate(x).map_err(to_py_err)
    }

    /// Get the degrees (n, m).
    fn degrees(&self) -> (usize, usize) {
        self.inner.degrees()
    }

    /// Get the number of subintervals.
    fn num_intervals(&self) -> usize {
        self.inner.num_intervals()
    }

    fn __repr__(&self) -> String {
        let (n, m) = self.inner.degrees();
        format!("PiecewiseRational([{}/{}], intervals={})", n, m, self.inner.num_intervals())
    }
}

/// Hermite interpolation constraints for rational approximants.
///
/// Exposes the constraint equations used in two-point Padé construction,
/// enabling custom solvers, optimization, and root-finding.
///
/// Parameters
/// ----------
/// left_derivs : array_like
///     Derivative values at left endpoint [f(a), f'(a), ...]
/// right_derivs : array_like
///     Derivative values at right endpoint [f(b), f'(b), ...]
/// n : int
///     Numerator degree
/// m : int
///     Denominator degree
/// a : float
///     Left endpoint
/// b : float
///     Right endpoint
///
/// Examples
/// --------
/// >>> import gelfgren as gf
/// >>> import numpy as np
/// >>> from scipy.optimize import least_squares
/// >>>
/// >>> # Create constraints for e^x on [0,1] with [2/1] rational
/// >>> left = np.array([1.0, 1.0])
/// >>> right = np.array([np.e, np.e])
/// >>> constraints = gf.HermiteConstraints(left, right, 2, 1, 0.0, 1.0)
/// >>>
/// >>> # Get linear system
/// >>> a, b, n = constraints.linear_system()
/// >>>
/// >>> # Or use for root-finding/optimization
/// >>> result = least_squares(constraints.residuals, x0=np.ones(4))
#[pyclass(name = "HermiteConstraints")]
struct PyHermiteConstraints {
    inner: CoreHermiteConstraints<f64>,
}

#[pymethods]
impl PyHermiteConstraints {
    #[new]
    fn new(
        left_derivs: PyReadonlyArray1<f64>,
        right_derivs: PyReadonlyArray1<f64>,
        n: usize,
        m: usize,
        a: f64,
        b: f64,
    ) -> PyResult<Self> {
        let left_vec = left_derivs.as_slice()?.to_vec();
        let right_vec = right_derivs.as_slice()?.to_vec();

        let inner = CoreHermiteConstraints::new(
            left_vec, right_vec, n, m, a, b
        ).map_err(to_py_err)?;

        Ok(Self { inner })
    }

    /// Get the linear system Ax = b representing the matching conditions.
    ///
    /// Returns
    /// -------
    /// A : ndarray, shape (2p, n+m+1)
    ///     Constraint matrix
    /// b : ndarray, shape (2p,)
    ///     Right-hand side vector
    /// n : int
    ///     Number of unknowns
    fn linear_system<'py>(&self, py: Python<'py>) -> (
        Bound<'py, PyArray1<f64>>,
        Bound<'py, PyArray1<f64>>,
        usize
    ) {
        let (a, b, n) = self.inner.linear_system();
        (
            a.to_pyarray_bound(py),
            b.to_pyarray_bound(py),
            n
        )
    }

    /// Evaluate constraint residuals given coefficient vector.
    ///
    /// Parameters
    /// ----------
    /// coeffs : array_like, shape (n+m+1,)
    ///     Coefficient vector [a₀, ..., aₙ, b₁, ..., bₘ]
    ///
    /// Returns
    /// -------
    /// residuals : ndarray, shape (2p,)
    ///     Residual vector r where rᵢ = R^(k)(xⱼ) - f^(k)(xⱼ)
    fn residuals<'py>(&self, py: Python<'py>, coeffs: PyReadonlyArray1<f64>) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let coeffs_vec = coeffs.as_slice()?.to_vec();
        let residuals = self.inner.residuals(&coeffs_vec).map_err(to_py_err)?;
        Ok(residuals.to_pyarray_bound(py))
    }

    /// Objective function: sum of squared residuals.
    ///
    /// Returns ||r||² where r is the residual vector.
    ///
    /// Parameters
    /// ----------
    /// coeffs : array_like, shape (n+m+1,)
    ///     Coefficient vector
    ///
    /// Returns
    /// -------
    /// float
    ///     Sum of squared residuals
    fn objective(&self, coeffs: PyReadonlyArray1<f64>) -> PyResult<f64> {
        let coeffs_vec = coeffs.as_slice()?.to_vec();
        self.inner.objective(&coeffs_vec).map_err(to_py_err)
    }

    /// Compute the Jacobian matrix of residuals.
    ///
    /// Parameters
    /// ----------
    /// coeffs : array_like, shape (n+m+1,)
    ///     Coefficient vector
    ///
    /// Returns
    /// -------
    /// jacobian : ndarray, shape (2p, n+m+1)
    ///     Jacobian matrix where J[i,j] = ∂rᵢ/∂cⱼ
    fn jacobian<'py>(&self, py: Python<'py>, coeffs: PyReadonlyArray1<f64>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let coeffs_vec = coeffs.as_slice()?.to_vec();
        let jac_flat = self.inner.jacobian(&coeffs_vec).map_err(to_py_err)?;
        let n_constraints = self.inner.num_constraints();
        let n_unknowns = self.inner.num_unknowns();

        // Reshape to 2D
        let jac_2d: Vec<Vec<f64>> = jac_flat
            .chunks(n_unknowns)
            .map(|row| row.to_vec())
            .collect();

        Ok(PyArray2::from_vec2_bound(py, &jac_2d)?)
    }

    /// Get quadratic form matrices for least-squares optimization.
    ///
    /// For minimizing ||Ax - b||², this is equivalent to minimizing
    /// x^T H x + g^T x where H = A^T A and g = -2 A^T b.
    ///
    /// Returns
    /// -------
    /// H : ndarray, shape (n+m+1, n+m+1)
    ///     Hessian matrix
    /// g : ndarray, shape (n+m+1,)
    ///     Gradient vector
    fn quadratic_form<'py>(&self, py: Python<'py>) -> (
        Bound<'py, PyArray2<f64>>,
        Bound<'py, PyArray1<f64>>
    ) {
        let (h_flat, g) = self.inner.quadratic_form();
        let n = self.inner.num_unknowns();

        // Reshape Hessian to 2D
        let h_2d: Vec<Vec<f64>> = h_flat
            .chunks(n)
            .map(|row| row.to_vec())
            .collect();

        (
            PyArray2::from_vec2_bound(py, &h_2d).unwrap(),
            g.to_pyarray_bound(py)
        )
    }

    /// Build rational function from coefficient vector.
    ///
    /// Parameters
    /// ----------
    /// coeffs : array_like, shape (n+m+1,)
    ///     Coefficient vector [a₀, ..., aₙ, b₁, ..., bₘ]
    ///
    /// Returns
    /// -------
    /// RationalFunction
    ///     Rational function R(x) = P(x)/Q(x)
    fn build_rational(&self, coeffs: PyReadonlyArray1<f64>) -> PyResult<PyRational> {
        let coeffs_vec = coeffs.as_slice()?.to_vec();
        let rational = self.inner.build_rational(&coeffs_vec).map_err(to_py_err)?;
        Ok(PyRational { inner: rational })
    }

    /// Number of constraint equations (2p).
    fn num_constraints(&self) -> usize {
        self.inner.num_constraints()
    }

    /// Number of unknowns (n+m+1).
    fn num_unknowns(&self) -> usize {
        self.inner.num_unknowns()
    }

    /// Numerator degree.
    fn numerator_degree(&self) -> usize {
        self.inner.numerator_degree()
    }

    /// Denominator degree.
    fn denominator_degree(&self) -> usize {
        self.inner.denominator_degree()
    }

    /// Order (number of derivatives at each endpoint).
    fn order(&self) -> usize {
        self.inner.order()
    }

    /// Interval [a, b].
    fn interval(&self) -> (f64, f64) {
        self.inner.interval()
    }

    fn __repr__(&self) -> String {
        let (n, m) = (self.inner.numerator_degree(), self.inner.denominator_degree());
        let p = self.inner.order();
        format!("HermiteConstraints([{}/{}], order={}, constraints={})", n, m, p, self.inner.num_constraints())
    }
}

/// Gelfgren: Piecewise Rational Interpolation Library
///
/// A high-performance library for piecewise rational interpolation and approximation,
/// implementing two-point Padé approximants.
///
/// Classes
/// -------
/// BernsteinPolynomial
///     Numerically stable polynomial representation
/// RationalFunction
///     Rational function P(x)/Q(x)
/// PadeApproximant
///     Padé approximants from power series
/// TwoPointPade
///     Two-point Padé approximants matching endpoint derivatives
/// PiecewiseRational
///     Piecewise rational approximation over mesh
/// Mesh
///     Mesh partitions for piecewise approximation
#[pymodule]
fn gelfgren(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyBernstein>()?;
    m.add_class::<PyRational>()?;
    m.add_class::<PyPade>()?;
    m.add_class::<PyTwoPointPade>()?;
    m.add_class::<PyPiecewiseRational>()?;
    m.add_class::<PyMesh>()?;
    m.add_class::<PyHermiteConstraints>()?;

    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add("__doc__", "Gelfgren: Piecewise Rational Interpolation Library")?;

    Ok(())
}
