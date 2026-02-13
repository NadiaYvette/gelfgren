//! Python bindings for Gelfgren library using PyO3.
//!
//! Provides Pythonic interface with NumPy integration and exception handling.

use gelfgren_core::bernstein::BernsteinPolynomial as CoreBernstein;
use gelfgren_core::rational::RationalFunction as CoreRational;
use gelfgren_core::pade::PadeApproximant as CorePade;
use gelfgren_core::piecewise::Mesh as CoreMesh;
use gelfgren_core::error::GelfgrenError;

use pyo3::prelude::*;
use pyo3::exceptions::{PyValueError, PyRuntimeError};
use numpy::{PyArray1, PyReadonlyArray1, ToPyArray};

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

    fn __repr__(&self) -> String {
        let (a, b) = self.inner.interval();
        format!("Mesh(interval=[{}, {}], points={})", a, b, self.inner.num_points())
    }
}

/// Gelfgren: Piecewise Rational Interpolation Library
///
/// A high-performance library for piecewise rational interpolation and approximation,
/// implementing algorithms from Jan Gelfgren's 1975 paper.
///
/// Classes
/// -------
/// BernsteinPolynomial
///     Numerically stable polynomial representation
/// RationalFunction
///     Rational function P(x)/Q(x)
/// PadeApproximant
///     Padé approximants from power series
/// Mesh
///     Mesh partitions for piecewise approximation
#[pymodule]
fn gelfgren(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyBernstein>()?;
    m.add_class::<PyRational>()?;
    m.add_class::<PyPade>()?;
    m.add_class::<PyMesh>()?;

    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add("__doc__", "Gelfgren: Piecewise Rational Interpolation Library")?;

    Ok(())
}
