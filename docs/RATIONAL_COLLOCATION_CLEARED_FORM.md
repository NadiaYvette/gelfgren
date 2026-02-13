# Rational Collocation with Cleared Denominators

## The Key Insight

**Question:** What if we multiply the differential equation by Q(x) to eliminate quotients before differentiating?

This transforms the problem and potentially provides significant advantages!

---

## Standard Formulation (Quotient Form)

### Approach

Given ODE: `-u'' = f(x)` with `u = P(x)/Q(x)`

**Directly differentiate the quotient:**
```
u' = (P'Q - PQ')/Q²
u'' = (P''Q² - 2P'Q'Q - PQ''Q + 2PQ'²)/Q³
```

**Collocation equation at xᵢ:**
```
-(P''Q² - 2P'Q'Q - PQ''Q + 2PQ'²)/Q³|ₓᵢ = f(xᵢ)
```

**Problems:**
- Complex quotient with Q³ in denominator
- Potential numerical instability if Q(xᵢ) is small
- Not obvious how to prevent spurious poles
- Difficult to analyze

---

## Cleared Formulation (Product Form)

### Derivation

**Start with:** `u = P/Q`

**Fundamental identity:** `Q·u = P`

**Key idea:** Multiply the ODE by Q³ to clear all denominators!

### Step-by-Step

#### 1. Express derivatives using product rule

From `Q·u = P`:
```
Differentiate once:
Q'·u + Q·u' = P'

Differentiate twice:
(Q'·u + Q·u')' = P''
Q''·u + Q'·u' + Q'·u' + Q·u'' = P''
Q''·u + 2Q'·u' + Q·u'' = P''
```

Solve for `Q·u''`:
```
Q·u'' = P'' - Q''·u - 2Q'·u'
```

#### 2. Substitute back into ODE

ODE: `-u'' = f`

Multiply by Q:
```
-Q·u'' = Q·f
```

Substitute expression for Q·u'':
```
-(P'' - Q''·u - 2Q'·u') = Q·f
-P'' + Q''·u + 2Q'·u' = Q·f
```

#### 3. Eliminate remaining u and u'

From `Q·u = P`: `u = P/Q`

From `Q'·u + Q·u' = P'`: `u' = (P' - Q'·u)/Q = (P' - Q'P/Q)/Q = (P'Q - Q'P)/Q²`

**Substitute into ODE:**
```
-P'' + Q''·(P/Q) + 2Q'·((P'Q - Q'P)/Q²) = Q·f
```

Multiply through by Q²:
```
-Q²·P'' + Q·Q''·P + 2Q'·(P'Q - Q'P) = Q³·f
-Q²·P'' + Q·Q''·P + 2Q·Q'·P' - 2Q'²·P = Q³·f
```

**Rearrange:**
```
Q²·P'' - Q·Q''·P - 2Q·Q'·P' + 2Q'²·P = -Q³·f
```

Or, collecting P terms:
```
Q²·P'' - 2Q·Q'·P' + (2Q'² - Q·Q'')·P = -Q³·f
```

---

## The Cleared Collocation Equation

### Final Form

At collocation point xᵢ:
```
Q²(xᵢ)·P''(xᵢ) - 2Q(xᵢ)·Q'(xᵢ)·P'(xᵢ) + [2Q'²(xᵢ) - Q(xᵢ)·Q''(xᵢ)]·P(xᵢ) = -Q³(xᵢ)·f(xᵢ)
```

**This is the cleared form!**

### Key Properties

**1. No quotients!**
- All terms are polynomials
- No Q in denominator
- More numerically stable

**2. Still nonlinear**
- Products of P and Q coefficients appear
- Terms like Q·Q', Q'², Q·Q'', Q³
- Still requires iterative solver

**3. Different structure**
- Linear operator acts on P: `Q²·P'' - 2Q·Q'·P' + (...)·P`
- Modulation by Q and its derivatives
- RHS is `-Q³·f`

**4. Spurious pole control**
- Equation involves Q³(xᵢ)
- If Q(xᵢ) → 0, the RHS vanishes proportionally
- Natural penalty against poles at collocation points!

---

## Comparison of Formulations

| Aspect | Quotient Form | Cleared Form |
|--------|---------------|--------------|
| **Derivatives** | u = P/Q, differentiate | Q·u = P, differentiate |
| **Final Equation** | Complex quotient with Q³ | Polynomial products only |
| **Numerical Stability** | Issues if Q(xᵢ) small | Better (no division) |
| **Pole Prevention** | Requires explicit constraints | Natural penalty via Q³ |
| **Nonlinearity** | Products + quotients | Products only |
| **Implementation** | Compute u''(xᵢ) | Evaluate polynomial products |
| **Symmetry** | Not obvious | Clearer structure |

---

## Mathematical Interpretation

### What Does It Mean?

**Quotient form:**
```
-u'' = f  where u = P/Q
```
"The second derivative of the rational equals f"

**Cleared form:**
```
Q²·P'' - 2Q·Q'·P' + (2Q'² - Q·Q'')·P = -Q³·f
```
"The numerator polynomial P satisfies a modified differential equation weighted by Q"

### Physical Interpretation

Think of Q as a **weighting function**:
- Where Q is large: equation weights P'' more heavily
- Where Q is small: equation weights lower derivatives
- The Q³·f term on RHS: forcing is amplified where denominator is large

This makes intuitive sense! If u = P/Q and Q is large, then P must vary more rapidly to produce the same u.

---

## Implementation

### Standard Collocation (Quotient Form)

```python
def residuals_quotient(coeffs):
    a, b = unpack(coeffs)  # P coefficients, Q coefficients
    P = build_polynomial(a)
    Q = build_polynomial(b)

    res = []
    for xi in collocation_points:
        # Compute u''(xi) using quotient rule
        Px = P.deriv(1)(xi)
        Pxx = P.deriv(2)(xi)
        Qx = Q.deriv(1)(xi)
        Qxx = Q.deriv(2)(xi)
        Qval = Q(xi)

        uxx = (Pxx*Qval**2 - 2*Px*Qx*Qval - P(xi)*Qxx*Qval + 2*P(xi)*Qx**2) / Qval**3

        res.append(-uxx - f(xi))

    return res
```

**Issues:**
- Division by Q³
- Unstable if Q(xi) ≈ 0
- Doesn't naturally prevent poles

### Cleared Collocation (Product Form)

```python
def residuals_cleared(coeffs):
    a, b = unpack(coeffs)  # P coefficients, Q coefficients
    P = build_polynomial(a)
    Q = build_polynomial(b)

    res = []
    for xi in collocation_points:
        # Evaluate P, P', P'' and Q, Q', Q''
        Pval = P(xi)
        Px = P.deriv(1)(xi)
        Pxx = P.deriv(2)(xi)

        Qval = Q(xi)
        Qx = Q.deriv(1)(xi)
        Qxx = Q.deriv(2)(xi)

        # Cleared form: Q²·P'' - 2Q·Q'·P' + (2Q'² - Q·Q'')·P = -Q³·f
        lhs = (Qval**2 * Pxx
               - 2*Qval*Qx * Px
               + (2*Qx**2 - Qval*Qxx) * Pval)

        rhs = -Qval**3 * f(xi)

        res.append(lhs - rhs)

    return res
```

**Advantages:**
- No division!
- If Q(xi) small, both LHS and RHS scale appropriately
- Natural penalty against Q(xi) = 0

---

## Numerical Stability Analysis

### Quotient Form: Catastrophic Cancellation

If Q(xᵢ) ≈ 0:
```
uxx = (Pxx*Q² - 2Px*Qx*Q - P*Qxx*Q + 2P*Qx²)/Q³
    = (small - small - small + small) / (very small)
    = (catastrophic cancellation) / (near zero)
```

**Result:** Complete loss of precision or NaN

### Cleared Form: Graceful Degradation

If Q(xᵢ) ≈ 0:
```
LHS = Q²·Pxx - 2Q·Qx·Px + (2Qx² - Q·Qxx)·P
    ≈ 2Qx²·P  (only Qx² term survives)

RHS = -Q³·f
    ≈ 0  (vanishes cubically)

Equation: 2Qx²·P ≈ 0
```

**Result:** Requires P(xᵢ) ≈ 0 if Qx ≠ 0
- **Natural pole prevention!**
- If P and Q both vanish: L'Hôpital can be used
- No division by zero

---

## Pole Prevention Mechanism

### Explicit Analysis

At collocation point xᵢ, the cleared equation is:
```
Q²(xᵢ)·P''(xᵢ) - 2Q(xᵢ)·Q'(xᵢ)·P'(xᵢ) + [2Q'²(xᵢ) - Q(xᵢ)·Q''(xᵢ)]·P(xᵢ) = -Q³(xᵢ)·f(xᵢ)
```

**Suppose Q(xᵢ) → 0:**

Leading order terms:
```
2Q'²(xᵢ)·P(xᵢ) ≈ 0  (since other terms have Q or Q² factors)
```

If Q'(xᵢ) ≠ 0, this forces:
```
P(xᵢ) ≈ 0
```

**Interpretation:** If denominator vanishes at xᵢ, numerator must also vanish!

**L'Hôpital applicable:**
```
u(xᵢ) = P(xᵢ)/Q(xᵢ) = 0/0  →  can evaluate as P'(xᵢ)/Q'(xᵢ)
```

This is a **removable singularity**, not a pole!

### Constraint Interpretation

The cleared form implicitly enforces:
```
If Q(xᵢ) = 0, then P(xᵢ) = 0  (assuming Q'(xᵢ) ≠ 0)
```

This is exactly the **compatibility condition** for having a finite limit!

No need for explicit constraints like "Q(x) > δ" - it's **built into the formulation**.

---

## Advantages Summary

### 1. **Numerical Stability**
- No division operations
- No catastrophic cancellation
- All terms are polynomial evaluations
- Standard floating-point arithmetic suffices

### 2. **Natural Pole Prevention**
- Spurious poles automatically penalized
- Equation enforces compatibility: Q=0 ⟹ P=0
- No need for inequality constraints
- Cleaner optimization landscape

### 3. **Clearer Structure**
- Polynomial products only
- Can be written in matrix form more easily
- Better conditioning possible
- Easier to analyze convergence

### 4. **Implementation Simplicity**
- No special case handling for small Q
- More robust to poor initial guesses
- Fewer numerical issues in practice

### 5. **Physical Interpretation**
- Q acts as a weighting function
- Natural variational interpretation
- Connection to weighted residual methods

---

## Potential Disadvantages

### 1. **Higher Polynomial Degrees**
Products like Q²·P'', Q³·f involve higher degree polynomials than P/Q itself.

**Impact:** May need more careful evaluation (but no fundamental issue)

### 2. **Still Nonlinear**
Products Q·Q', Q'², Q·Q'', Q³ still appear.

**Impact:** Still requires iterative solver (but quotient form does too!)

### 3. **Boundary Conditions Need Care**
BCs u(a) = α becomes Q(a)·α = P(a), which is linear in both P and Q.

**Impact:** Need to normalize (e.g., fix Q(a) = 1) or handle as constraint

### 4. **Less Familiar**
Standard collocation literature uses quotient form.

**Impact:** May need to develop new theory for convergence analysis

---

## Practical Example

### Problem
```
-u'' = 2 on [0,1]
u(0) = 0, u(1) = 0
```

Exact solution: `u(x) = x(1-x)`

### Using [2/1] Rational: u = (a₀ + a₁x + a₂x²)/(1 + b₁x)

#### Quotient Form

```python
def residuals_quotient(coeffs):
    a0, a1, a2, b1 = coeffs
    res = []

    # BC: u(0) = a0/1 = 0
    res.append(a0)

    # BC: u(1) = (a0+a1+a2)/(1+b1) = 0
    res.append(a0 + a1 + a2)  # numerator = 0

    # Collocation at x=1/3, x=2/3
    for x in [1/3, 2/3]:
        P = a0 + a1*x + a2*x**2
        Px = a1 + 2*a2*x
        Pxx = 2*a2
        Q = 1 + b1*x
        Qx = b1
        Qxx = 0

        uxx = (Pxx*Q**2 - 2*Px*Qx*Q - P*Qxx*Q + 2*P*Qx**2) / Q**3
        res.append(-uxx - 2.0)

    return res
```

#### Cleared Form

```python
def residuals_cleared(coeffs):
    a0, a1, a2, b1 = coeffs
    res = []

    # BC: Q(0)·u(0) = P(0)
    # 1·0 = a0  →  a0 = 0
    res.append(a0)

    # BC: Q(1)·u(1) = P(1)
    # (1+b1)·0 = a0+a1+a2  →  a0+a1+a2 = 0
    res.append(a0 + a1 + a2)

    # Collocation at x=1/3, x=2/3
    for x in [1/3, 2/3]:
        P = a0 + a1*x + a2*x**2
        Px = a1 + 2*a2*x
        Pxx = 2*a2
        Q = 1 + b1*x
        Qx = b1
        Qxx = 0

        # Q²·P'' - 2Q·Qx·Px + (2Qx² - Q·Qxx)·P = -Q³·f
        lhs = Q**2 * Pxx - 2*Q*Qx * Px + (2*Qx**2 - Q*Qxx) * P
        rhs = -Q**3 * 2.0

        res.append(lhs - rhs)

    return res
```

### Running Both

```python
from scipy.optimize import fsolve

# Quotient form
x0 = [0, 1, -1, 0]
sol_q = fsolve(residuals_quotient, x0)
print(f"Quotient: {sol_q}")

# Cleared form
sol_c = fsolve(residuals_cleared, x0)
print(f"Cleared: {sol_c}")

# Should both converge to same solution
# Expected: a0=0, a1=1, a2=-1, b1≈0 (reduces to polynomial x - x²)
```

**Expected behavior:**
- Both should find same solution
- Cleared form may converge faster
- Cleared form more robust to poor initial guess

---

## Connection to Weighted Residual Methods

### Standard Weighted Residuals

General form:
```
∫ w(x)·[L[u] - f]·ψ(x) dx = 0
```

where w(x) is a weight function.

### Cleared Collocation as Weighted Residual

Our cleared form:
```
Q²·P'' - 2Q·Qx·Px + (2Qx² - Q·Qxx)·P = -Q³·f
```

can be written as:
```
Q²·[L[P]] = -Q³·f
```

where L is a differential operator.

Rearranging:
```
L[P] + Q·f = 0
```

This is a **weighted residual** with weight function Q!

**Interpretation:** We're solving the residual equation with denominator as weight.

This provides a variational framework for analyzing the method.

---

## Theoretical Advantages for Convergence Analysis

### 1. **Operator Structure**

Cleared form has structure:
```
L_Q[P] = -Q³·f
```

where L_Q is a differential operator depending on Q.

This is closer to standard PDE theory than quotient form.

### 2. **Energy Estimates**

Can potentially multiply by P and integrate:
```
∫ P·L_Q[P] dx = ∫ P·(-Q³·f) dx
```

Left side may have coercivity properties useful for existence/uniqueness proofs.

### 3. **Fixed Point Formulation**

Can solve iteratively:
1. Fix Q, solve for P (linear problem!)
2. Fix P, update Q (nonlinear, but one variable)
3. Iterate

**Quotient form doesn't separate as cleanly.**

---

## Recommendations

### When to Use Cleared Form

**Strongly recommended:**
1. When Q might have small values
2. Near-pole problems where stability is critical
3. Optimization-based solvers (smoother landscape)
4. Need for robust implementation

**Consider cleared form:**
1. When developing new rational collocation software
2. For production code requiring reliability
3. When theory development is needed

### When Quotient Form Acceptable

1. Q guaranteed bounded away from zero
2. Only comparing with existing literature
3. Quick prototype implementation
4. When Q = 1 + ε·Q̃ with small ε (perturbative)

---

## Open Questions

### Research Directions

1. **Convergence theory:** Does cleared form have better convergence properties?

2. **Optimal Q strategy:** How to choose Q for best conditioning?

3. **Preconditioning:** Can cleared form lead to better preconditioners?

4. **Extensions:** How does this generalize to:
   - Higher-order equations?
   - Systems of ODEs?
   - PDEs?

5. **Comparison:** Numerical experiments on standard test problems

6. **Hybrid methods:** Can we combine cleared and quotient forms adaptively?

---

## Implementation in Gelfgren

### Proposed API

```rust
pub enum CollocationFormulation {
    Quotient,   // Standard u = P/Q form
    Cleared,    // Q²·P'' form (recommended)
}

impl RationalCollocation {
    pub fn new(formulation: CollocationFormulation) -> Self { ... }

    pub fn solve(&self, ode: ODE, boundary_conditions: BC) -> RationalFunction {
        match self.formulation {
            Quotient => self.solve_quotient_form(ode, boundary_conditions),
            Cleared => self.solve_cleared_form(ode, boundary_conditions),
        }
    }
}
```

### Python Example

```python
import gelfgren as gf

# Create collocation solver with cleared formulation
solver = gf.RationalCollocation(
    formulation='cleared',  # or 'quotient'
    n=2, m=1,  # [2/1] rational
    collocation_points=[1/3, 2/3]
)

# Solve BVP
solution = solver.solve(
    ode=lambda x: 2.0,  # -u'' = 2
    boundary_conditions=((0, 0.0), (1, 0.0))
)

print(solution.evaluate(0.5))  # Should be close to 0.25
```

---

## Conclusion

**Key Insight:** Clearing the denominator before differentiating transforms rational collocation from:
- Unstable quotient with potential division by zero
- To stable polynomial products with built-in pole prevention

**Main advantages:**
1. ✅ Numerical stability (no division)
2. ✅ Natural pole prevention (compatibility enforced)
3. ✅ Clearer mathematical structure
4. ✅ Better for optimization algorithms

**Recommendation:** **Use cleared form for production implementation!**

The cleared formulation is essentially a **weighted residual method** where the weight is the denominator polynomial. This provides both practical benefits (stability) and theoretical advantages (variational framework).

---

## References

1. **Brezinski, C., Redivo-Zaglia, M.** (1991). *Extrapolation Methods: Theory and Practice*. North-Holland.
   - Discusses clearing denominators in rational approximation

2. **Saff, E.B., Totik, V.** (1997). *Logarithmic Potentials with External Fields*. Springer.
   - Theoretical foundations for rational approximation

3. **Boyd, J.P.** (2001). *Chebyshev and Fourier Spectral Methods*. Dover.
   - Chapter 17: Rational Chebyshev functions
   - Discusses pole management strategies

4. **Berrut, J.P., Mittelmann, H.D.** (2004). "Matrices for the direct determination of the barycentric weights of rational interpolation." *J. Comput. Appl. Math.* 164-165: 19-28.
   - Stability considerations in rational approximation

---

**Document Status:** ✅ Complete analysis of cleared denominator formulation

**Related Documents:**
- `RATIONAL_COLLOCATION_METHODS.md` - General collocation theory
- `CONSTRAINT_INTERFACE_COMPLETE.md` - HermiteConstraints implementation
- `BENCHMARK_REPORT_UPDATE.md` - BVP methodology clarification

**Next Steps:** Implement both formulations in Gelfgren and benchmark stability/accuracy!
