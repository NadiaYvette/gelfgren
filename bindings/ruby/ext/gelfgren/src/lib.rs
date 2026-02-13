use gelfgren_core::{
    bernstein::BernsteinPolynomial as CoreBernstein,
    rational::RationalFunction as CoreRational,
    pade::PadeApproximant as CorePade,
    error::GelfgrenError,
};
use magnus::{
    class, define_class, define_module, function, method, prelude::*, Error, Module, Object,
    RArray, RClass, Ruby, TryConvert, Value,
};

/// Wrapper for BernsteinPolynomial
#[derive(Clone)]
#[magnus::wrap(class = "Gelfgren::BernsteinPolynomial")]
struct BernsteinPolynomial {
    inner: CoreBernstein<f64>,
}

impl BernsteinPolynomial {
    fn new(coeffs: RArray, a: f64, b: f64) -> Result<Self, Error> {
        let coeffs_vec: Vec<f64> = coeffs
            .each()
            .map(|v| f64::try_convert(v?))
            .collect::<Result<Vec<_>, _>>()?;

        CoreBernstein::new(coeffs_vec, a, b)
            .map(|inner| Self { inner })
            .map_err(|e| Error::new(magnus::exception::runtime_error(), format!("{:?}", e)))
    }

    fn evaluate(&self, x: Value) -> Result<Value, Error> {
        let ruby = Ruby::get_with(x);

        // Check if x is an array
        if let Ok(arr) = RArray::try_convert(x) {
            // Vectorized evaluation
            let results: Result<Vec<_>, _> = arr
                .each()
                .map(|v| {
                    let xi: f64 = v?;
                    Ok(self.inner.evaluate(xi))
                })
                .collect();

            let results = results?;
            Ok(RArray::from_vec(results).as_value())
        } else {
            // Single value evaluation
            let xi: f64 = x.try_convert()?;
            let result = self.inner.evaluate(xi);
            Ok(ruby.into_value(result))
        }
    }

    fn derivative(&self) -> Result<Self, Error> {
        Ok(Self {
            inner: self.inner.derivative(),
        })
    }

    fn integral(&self) -> Result<Self, Error> {
        Ok(Self {
            inner: self.inner.integral(),
        })
    }

    fn degree(&self) -> usize {
        self.inner.degree()
    }

    fn interval(&self) -> RArray {
        RArray::from_vec(vec![self.inner.a(), self.inner.b()])
    }

    fn add(&self, other: &Self) -> Result<Self, Error> {
        (self.inner.clone() + other.inner.clone())
            .map(|inner| Self { inner })
            .map_err(|e| Error::new(magnus::exception::runtime_error(), format!("{:?}", e)))
    }

    fn subtract(&self, other: &Self) -> Result<Self, Error> {
        (self.inner.clone() - other.inner.clone())
            .map(|inner| Self { inner })
            .map_err(|e| Error::new(magnus::exception::runtime_error(), format!("{:?}", e)))
    }

    fn multiply(&self, other: &Self) -> Result<Self, Error> {
        (self.inner.clone() * other.inner.clone())
            .map(|inner| Self { inner })
            .map_err(|e| Error::new(magnus::exception::runtime_error(), format!("{:?}", e)))
    }

    fn scale(&self, scalar: f64) -> Self {
        Self {
            inner: self.inner.clone() * scalar,
        }
    }

    fn elevate_degree(&self, r: usize) -> Self {
        Self {
            inner: self.inner.elevate_degree(r),
        }
    }

    fn to_s(&self) -> String {
        format!(
            "BernsteinPolynomial(degree={}, interval=[{}, {}])",
            self.inner.degree(),
            self.inner.a(),
            self.inner.b()
        )
    }
}

/// Wrapper for RationalFunction
#[derive(Clone)]
#[magnus::wrap(class = "Gelfgren::RationalFunction")]
struct RationalFunction {
    inner: CoreRational<f64>,
}

impl RationalFunction {
    fn new(numerator: &BernsteinPolynomial, denominator: &BernsteinPolynomial) -> Result<Self, Error> {
        CoreRational::new(numerator.inner.clone(), denominator.inner.clone())
            .map(|inner| Self { inner })
            .map_err(|e| Error::new(magnus::exception::runtime_error(), format!("{:?}", e)))
    }

    fn evaluate(&self, x: Value) -> Result<Value, Error> {
        let ruby = Ruby::get_with(x);

        if let Ok(arr) = RArray::try_convert(x) {
            // Vectorized evaluation
            let results: Result<Vec<_>, _> = arr
                .each()
                .map(|v| {
                    let xi: f64 = v?;
                    self.inner
                        .evaluate(xi)
                        .map_err(|e| Error::new(magnus::exception::runtime_error(), format!("{:?}", e)))
                })
                .collect();

            let results = results?;
            Ok(RArray::from_vec(results).as_value())
        } else {
            // Single value evaluation
            let xi: f64 = x.try_convert()?;
            self.inner
                .evaluate(xi)
                .map(|result| ruby.into_value(result))
                .map_err(|e| Error::new(magnus::exception::runtime_error(), format!("{:?}", e)))
        }
    }

    fn derivative(&self) -> Self {
        Self {
            inner: self.inner.derivative(),
        }
    }

    fn numerator(&self) -> BernsteinPolynomial {
        BernsteinPolynomial {
            inner: self.inner.numerator().clone(),
        }
    }

    fn denominator(&self) -> BernsteinPolynomial {
        BernsteinPolynomial {
            inner: self.inner.denominator().clone(),
        }
    }

    fn add(&self, other: &Self) -> Result<Self, Error> {
        (self.inner.clone() + other.inner.clone())
            .map(|inner| Self { inner })
            .map_err(|e| Error::new(magnus::exception::runtime_error(), format!("{:?}", e)))
    }

    fn subtract(&self, other: &Self) -> Result<Self, Error> {
        (self.inner.clone() - other.inner.clone())
            .map(|inner| Self { inner })
            .map_err(|e| Error::new(magnus::exception::runtime_error(), format!("{:?}", e)))
    }

    fn multiply(&self, other: &Self) -> Result<Self, Error> {
        (self.inner.clone() * other.inner.clone())
            .map(|inner| Self { inner })
            .map_err(|e| Error::new(magnus::exception::runtime_error(), format!("{:?}", e)))
    }

    fn divide(&self, other: &Self) -> Result<Self, Error> {
        (self.inner.clone() / other.inner.clone())
            .map(|inner| Self { inner })
            .map_err(|e| Error::new(magnus::exception::runtime_error(), format!("{:?}", e)))
    }

    fn reciprocal(&self) -> Self {
        Self {
            inner: self.inner.reciprocal(),
        }
    }

    fn to_s(&self) -> String {
        format!(
            "RationalFunction(num_degree={}, den_degree={}, interval=[{}, {}])",
            self.inner.numerator().degree(),
            self.inner.denominator().degree(),
            self.inner.numerator().a(),
            self.inner.numerator().b()
        )
    }
}

/// Wrapper for PadeApproximant
#[derive(Clone)]
#[magnus::wrap(class = "Gelfgren::PadeApproximant")]
struct PadeApproximant {
    inner: CorePade<f64>,
}

impl PadeApproximant {
    fn from_power_series(coeffs: RArray, n: usize, m: usize, a: f64, b: f64) -> Result<Self, Error> {
        let coeffs_vec: Vec<f64> = coeffs
            .each()
            .map(|v| f64::try_convert(v?))
            .collect::<Result<Vec<_>, _>>()?;

        CorePade::from_power_series(coeffs_vec, n, m, a, b)
            .map(|inner| Self { inner })
            .map_err(|e| Error::new(magnus::exception::runtime_error(), format!("{:?}", e)))
    }

    fn from_derivatives(derivs: RArray, n: usize, m: usize, x0: f64, a: f64, b: f64) -> Result<Self, Error> {
        let derivs_vec: Vec<f64> = derivs
            .each()
            .map(|v| f64::try_convert(v?))
            .collect::<Result<Vec<_>, _>>()?;

        CorePade::from_derivatives(derivs_vec, n, m, x0, a, b)
            .map(|inner| Self { inner })
            .map_err(|e| Error::new(magnus::exception::runtime_error(), format!("{:?}", e)))
    }

    fn evaluate(&self, x: Value) -> Result<Value, Error> {
        let ruby = Ruby::get_with(x);

        if let Ok(arr) = RArray::try_convert(x) {
            // Vectorized evaluation
            let results: Result<Vec<_>, _> = arr
                .each()
                .map(|v| {
                    let xi: f64 = v?;
                    self.inner
                        .evaluate(xi)
                        .map_err(|e| Error::new(magnus::exception::runtime_error(), format!("{:?}", e)))
                })
                .collect();

            let results = results?;
            Ok(RArray::from_vec(results).as_value())
        } else {
            // Single value evaluation
            let xi: f64 = x.try_convert()?;
            self.inner
                .evaluate(xi)
                .map(|result| ruby.into_value(result))
                .map_err(|e| Error::new(magnus::exception::runtime_error(), format!("{:?}", e)))
        }
    }

    fn rational(&self) -> RationalFunction {
        RationalFunction {
            inner: self.inner.rational().clone(),
        }
    }

    fn orders(&self) -> RArray {
        RArray::from_vec(vec![self.inner.n(), self.inner.m()])
    }

    fn to_s(&self) -> String {
        format!(
            "PadeApproximant([{}/{}], interval=[{}, {}])",
            self.inner.n(),
            self.inner.m(),
            self.inner.rational().numerator().a(),
            self.inner.rational().numerator().b()
        )
    }
}

#[magnus::init]
fn init(ruby: &Ruby) -> Result<(), Error> {
    // Create Gelfgren module
    let module = define_module("Gelfgren")?;

    // Define BernsteinPolynomial class
    let bernstein_class = module.define_class("BernsteinPolynomial", class::object())?;
    bernstein_class.define_singleton_method("new", function!(BernsteinPolynomial::new, 3))?;
    bernstein_class.define_method("evaluate", method!(BernsteinPolynomial::evaluate, 1))?;
    bernstein_class.define_method("derivative", method!(BernsteinPolynomial::derivative, 0))?;
    bernstein_class.define_method("integral", method!(BernsteinPolynomial::integral, 0))?;
    bernstein_class.define_method("degree", method!(BernsteinPolynomial::degree, 0))?;
    bernstein_class.define_method("interval", method!(BernsteinPolynomial::interval, 0))?;
    bernstein_class.define_method("+", method!(BernsteinPolynomial::add, 1))?;
    bernstein_class.define_method("-", method!(BernsteinPolynomial::subtract, 1))?;
    bernstein_class.define_method("*", method!(BernsteinPolynomial::multiply, 1))?;
    bernstein_class.define_method("scale", method!(BernsteinPolynomial::scale, 1))?;
    bernstein_class.define_method("elevate_degree", method!(BernsteinPolynomial::elevate_degree, 1))?;
    bernstein_class.define_method("to_s", method!(BernsteinPolynomial::to_s, 0))?;

    // Define RationalFunction class
    let rational_class = module.define_class("RationalFunction", class::object())?;
    rational_class.define_singleton_method("new", function!(RationalFunction::new, 2))?;
    rational_class.define_method("evaluate", method!(RationalFunction::evaluate, 1))?;
    rational_class.define_method("derivative", method!(RationalFunction::derivative, 0))?;
    rational_class.define_method("numerator", method!(RationalFunction::numerator, 0))?;
    rational_class.define_method("denominator", method!(RationalFunction::denominator, 0))?;
    rational_class.define_method("+", method!(RationalFunction::add, 1))?;
    rational_class.define_method("-", method!(RationalFunction::subtract, 1))?;
    rational_class.define_method("*", method!(RationalFunction::multiply, 1))?;
    rational_class.define_method("/", method!(RationalFunction::divide, 1))?;
    rational_class.define_method("reciprocal", method!(RationalFunction::reciprocal, 0))?;
    rational_class.define_method("to_s", method!(RationalFunction::to_s, 0))?;

    // Define PadeApproximant class
    let pade_class = module.define_class("PadeApproximant", class::object())?;
    pade_class.define_singleton_method("from_power_series", function!(PadeApproximant::from_power_series, 5))?;
    pade_class.define_singleton_method("from_derivatives", function!(PadeApproximant::from_derivatives, 6))?;
    pade_class.define_method("evaluate", method!(PadeApproximant::evaluate, 1))?;
    pade_class.define_method("rational", method!(PadeApproximant::rational, 0))?;
    pade_class.define_method("orders", method!(PadeApproximant::orders, 0))?;
    pade_class.define_method("to_s", method!(PadeApproximant::to_s, 0))?;

    Ok(())
}
