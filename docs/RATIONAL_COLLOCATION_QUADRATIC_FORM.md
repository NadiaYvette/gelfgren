# Quadratic Formulation: Rational Collocation with Explicit Solution Values

## The Brilliant Insight

**Question:** Instead of eliminating u, what if we treat u and u' as explicit unknowns at collocation points?

**Result:** The system becomes **purely quadratic** instead of cubic/higher-order!

---

## The Formulation

### Starting Point

Given ODE: `-u'' = f(x)`

Define: **P = Q·u** (P is the product of denominator and solution)

### Differentiate the Product

```
P = Q·u

P' = Q'·u + Q·u'  (product rule)

P'' = (Q'·u + Q·u')'
    = Q''·u + Q'·u' + Q'·u' + Q·u''
    = Q''·u + 2Q'·u' + Q·u''
```

### Substitute ODE

From `-u'' = f(x)`, we have `u'' = -f(x)`:

```
P'' = Q''·u + 2Q'·u' + Q·(-f(x))
P'' = Q''·u + 2Q'·u' - Q·f(x)
```

### The Three Coupled Equations

At each collocation point xᵢ:

```
(1)  P(xᵢ) = Q(xᵢ)·u(xᵢ)

(2)  P'(xᵢ) = Q'(xᵢ)·u(xᵢ) + Q(xᵢ)·u'(xᵢ)

(3)  P''(xᵢ) = Q''(xᵢ)·u(xᵢ) + 2Q'(xᵢ)·u'(xᵢ) - Q(xᵢ)·f(xᵢ)
```

**Key observation:** These are **linear** in P, P', P'', Q, Q', Q'', u, u' individually, but involve **products** of unknowns!

---

## The Quadratic Structure

### Unknowns

**Global unknowns (polynomial coefficients):**
- P coefficients: `{a₀, a₁, ..., aₙ}` → (n+1) unknowns
- Q coefficients: `{b₀, b₁, ..., bₘ}` → m unknowns (normalize b₀ = 1)

**Local unknowns (per collocation point xᵢ):**
- Solution value: `u(xᵢ)` → 1 unknown per point
- Solution derivative: `u'(xᵢ)` → 1 unknown per point

**Total for k collocation points:**
```
Total unknowns = (n+1) + m + 2k
```

### Equations

**Per collocation point xᵢ (3 equations):**
```
(1)  P(xᵢ) - Q(xᵢ)·u(xᵢ) = 0
(2)  P'(xᵢ) - Q'(xᵢ)·u(xᵢ) - Q(xᵢ)·u'(xᵢ) = 0
(3)  P''(xᵢ) - Q''(xᵢ)·u(xᵢ) - 2Q'(xᵢ)·u'(xᵢ) + Q(xᵢ)·f(xᵢ) = 0
```

**Plus boundary conditions (2 equations):**
```
(BC1)  u(a) = α
(BC2)  u(b) = β
```

**Total equations:**
```
Total = 3k + 2
```

### Counting Arguments

For well-determined system:
```
3k + 2 = (n+1) + m + 2k
k + 2 = n + m + 1
k = n + m - 1
```

**Example:** [2/1] rational (n=2, m=1)
```
k = 2 + 1 - 1 = 2 collocation points
Total unknowns: 3 + 1 + 4 = 8
Total equations: 3·2 + 2 = 8 ✓
```

---

## Why This Is QUADRATIC

### Analysis of Nonlinearity

Look at equation (1): `P(xᵢ) = Q(xᵢ)·u(xᵢ)`

**Left side:**
```
P(xᵢ) = a₀·B₀(xᵢ) + a₁·B₁(xᵢ) + ... + aₙ·Bₙ(xᵢ)
```
Linear in coefficients {a₀, ..., aₙ}

**Right side:**
```
Q(xᵢ)·u(xᵢ) = [b₁·B₁(xᵢ) + ... + bₘ·Bₘ(xᵢ)] · u(xᵢ)  (with b₀=1)
           = b₁·B₁(xᵢ)·u(xᵢ) + ... + bₘ·Bₘ(xᵢ)·u(xᵢ) + B₀(xᵢ)·u(xᵢ)
```

**Products of unknowns:**
- `b₁ × u(xᵢ)` - product of two unknowns → **quadratic term**
- `b₂ × u(xᵢ)` - product of two unknowns → **quadratic term**
- etc.

**Degree of polynomial system:** 2 (quadratic)

### Comparison with Previous Formulations

**Cleared form (eliminated u):**
```
Q²·P'' - 2Q·Q'·P' + (2Q'² - Q·Q'')·P = -Q³·f
```

Terms like:
- `Q·Q'·P'` - involves `bᵢ × bⱼ × aₖ` → cubic
- `Q'²·P` - involves `bᵢ × bⱼ × aₖ` → cubic
- `Q³·f` - involves `bᵢ × bⱼ × bₖ` → cubic

**Degree:** 3 (cubic)

**New form (explicit u):**
```
P - Q·u = 0
P' - Q'·u - Q·u' = 0
P'' - Q''·u - 2Q'·u' + Q·f = 0
```

All terms involve products of at most 2 unknowns:
- `Q·u` → quadratic
- `Q'·u` → quadratic
- `Q·u'` → quadratic

**Degree:** 2 (quadratic)

**QUADRATIC IS MUCH BETTER THAN CUBIC!**

---

## Advantages of Quadratic Structure

### 1. **Better Convergence Properties**

Quadratic systems have special structure:
- Newton's method converges quadratically (if close to solution)
- Better basin of attraction
- More predictable behavior

Cubic/higher-order systems:
- May have multiple spurious solutions
- Worse conditioning
- Harder to guarantee convergence

### 2. **Specialized Solvers Available**

Can use **quadratic programming** methods:
- Sequential quadratic programming (SQP)
- Trust-region methods designed for quadratic
- Levenberg-Marquardt (natural fit)

### 3. **Clearer Hessian Structure**

Second derivatives are simpler:
- Quadratic terms → Hessian has constant entries
- Cubic terms → Hessian depends on solution

**Better for optimization algorithms!**

### 4. **Natural Decoupling**

Can solve in stages:
1. Fix Q coefficients, solve for P and u values (linear!)
2. Update Q coefficients
3. Iterate

This is a **bilinear alternating** strategy.

### 5. **Physical Interpretation**

The unknowns have clear meaning:
- `{aᵢ}` - numerator shape
- `{bⱼ}` - denominator shape
- `u(xᵢ)` - actual solution values
- `u'(xᵢ)` - actual solution slopes

**Direct connection to physics!**

---

## Implementation

### Residual Function

```python
def residuals_quadratic(coeffs):
    """
    Quadratic formulation with explicit u values.

    coeffs = [a0, ..., an, b1, ..., bm, u1, u'1, u2, u'2, ..., uk, u'k]
    """
    # Unpack coefficients
    n_p = n + 1  # P coefficients
    n_q = m      # Q coefficients (b0 = 1 normalized)
    n_u = 2 * k  # u and u' at k collocation points

    a = coeffs[0:n_p]                    # P coefficients
    b = coeffs[n_p:n_p+n_q]              # Q coefficients (without b0)
    u_vals = coeffs[n_p+n_q::2]          # u values at collocation points
    up_vals = coeffs[n_p+n_q+1::2]       # u' values at collocation points

    # Build polynomials P and Q
    P = BernsteinPolynomial(a, x_left, x_right)
    Q = BernsteinPolynomial([1.0] + list(b), x_left, x_right)

    residuals = []

    # Boundary conditions
    residuals.append(u_vals[0] - alpha)   # u(a) = alpha
    residuals.append(u_vals[-1] - beta)   # u(b) = beta

    # Collocation equations
    for i, xi in enumerate(collocation_points):
        ui = u_vals[i]
        upi = up_vals[i]

        # Equation 1: P(xi) = Q(xi)·u(xi)
        res1 = P(xi) - Q(xi) * ui
        residuals.append(res1)

        # Equation 2: P'(xi) = Q'(xi)·u(xi) + Q(xi)·u'(xi)
        res2 = P.deriv(1)(xi) - Q.deriv(1)(xi) * ui - Q(xi) * upi
        residuals.append(res2)

        # Equation 3: P''(xi) = Q''(xi)·u(xi) + 2Q'(xi)·u'(xi) - Q(xi)·f(xi)
        res3 = (P.deriv(2)(xi)
                - Q.deriv(2)(xi) * ui
                - 2 * Q.deriv(1)(xi) * upi
                + Q(xi) * f(xi))
        residuals.append(res3)

    return np.array(residuals)
```

**Note:** All products are explicit `Q(xi) * ui`, not hidden in quotient rule!

### Solving

```python
from scipy.optimize import fsolve, least_squares

# Initial guess
x0 = np.zeros(n_p + n_q + 2*k)
# ... initialize with polynomial solution or interpolation ...

# Solve quadratic system
solution = fsolve(residuals_quadratic, x0)

# Or use least_squares for robustness
solution = least_squares(residuals_quadratic, x0, method='lm')
```

**Levenberg-Marquardt (method='lm')** is particularly good for quadratic systems!

---

## Practical Example

### Problem
```
-u'' = 2 on [0,1]
u(0) = 0, u(1) = 0
```

Exact: `u(x) = x(1-x)`

### Using [2/1] Rational

**Setup:**
- n = 2, m = 1 → k = 2 collocation points
- Choose x₁ = 1/3, x₂ = 2/3
- Unknowns: {a₀, a₁, a₂, b₁, u₁, u'₁, u₂, u'₂} → 8 unknowns
- Equations: 2 BC + 3×2 collocation → 8 equations ✓

**Implementation:**

```python
def residuals(coeffs):
    a0, a1, a2, b1, u1, up1, u2, up2 = coeffs

    res = []

    # BC: u(0) = 0 (boundary value at x=0, which is u1 if we include it)
    # Actually, for interior collocation, u(0) comes from P(0)/Q(0)
    # Let's say u1 is at x=1/3, u2 at x=2/3
    # We need to enforce BC differently...

    # At x=0: P(0) = Q(0)·0 → a0 = 0
    res.append(a0)

    # At x=1: P(1) = Q(1)·0 → a0+a1+a2 = 0
    res.append(a0 + a1 + a2)

    # Collocation at x1 = 1/3
    x = 1/3
    P = a0 + a1*x + a2*x**2
    Px = a1 + 2*a2*x
    Pxx = 2*a2
    Q = 1 + b1*x
    Qx = b1
    Qxx = 0

    # (1) P(x1) = Q(x1)·u1
    res.append(P - Q * u1)

    # (2) P'(x1) = Q'(x1)·u1 + Q(x1)·u'1
    res.append(Px - Qx * u1 - Q * up1)

    # (3) P''(x1) = Q''(x1)·u1 + 2Q'(x1)·u'1 - Q(x1)·f(x1)
    res.append(Pxx - Qxx * u1 - 2*Qx * up1 + Q * 2.0)

    # Collocation at x2 = 2/3
    x = 2/3
    P = a0 + a1*x + a2*x**2
    Px = a1 + 2*a2*x
    Pxx = 2*a2
    Q = 1 + b1*x
    Qx = b1
    Qxx = 0

    res.append(P - Q * u2)
    res.append(Px - Qx * u2 - Q * up2)
    res.append(Pxx - Qxx * u2 - 2*Qx * up2 + Q * 2.0)

    return res

# Solve
x0 = [0, 1, -1, 0, 0.22, 0.56, 0.22, -0.56]  # Initial guess
solution = fsolve(residuals, x0)

print(f"P coefficients: {solution[0:3]}")
print(f"Q coefficients: [1.0, {solution[3]}]")
print(f"u(1/3) = {solution[4]:.6f}, exact = {(1/3)*(2/3):.6f}")
print(f"u(2/3) = {solution[6]:.6f}, exact = {(2/3)*(1/3):.6f}")
```

---

## Comparison Table

| Formulation | Degree | Division? | Structure | Solver |
|-------------|--------|-----------|-----------|--------|
| **Quotient** | 3+ | Yes (/Q³) | Complex quotients | Newton, opt |
| **Cleared** | 3 | No | Polynomial products | Newton, opt |
| **Quadratic** | **2** | **No** | **Bilinear products** | **SQP, L-M** |

**Quadratic is best on all counts!**

---

## Theoretical Advantages

### 1. Bilinear Alternating Strategy

**Fix Q, solve for P and u:**
```
P(xi) = Q(xi)·u(xi)           [linear in P, u]
P'(xi) = Q'(xi)·u(xi) + Q(xi)·u'(xi)  [linear in P, u, u']
P''(xi) = Q''(xi)·u(xi) + 2Q'(xi)·u'(xi) - Q(xi)·f  [linear in P, u, u']
```

With Q fixed, this is a **linear system** in {a₀, ..., aₙ, u₁, u'₁, ..., uₖ, u'ₖ}!

**Fix P and u, solve for Q:**
```
P(xi) - Q(xi)·u(xi) = 0  [linear in Q if P, u fixed]
```

This enables **alternating optimization**:
1. Fix Q → solve for P, u (linear system)
2. Fix P, u → solve for Q (linear system)
3. Iterate until convergence

Each step is a **linear solve** - no Newton iteration needed!

### 2. Connection to Constrained Least Squares

Can reformulate as:
```
minimize ||P - Q·u||² over (P, Q, u)
subject to: ODE constraints
```

This is a **quadratically constrained quadratic program (QCQP)**!

Specialized solvers exist for QCQP.

### 3. Natural Regularization

Can add regularization naturally:
```
minimize ||P - Q·u||² + λ₁||Q||² + λ₂||u''||²
```

- λ₁: prevents Q from becoming too large
- λ₂: enforces smoothness of u

**This is much harder with eliminated-u formulations!**

---

## Spurious Pole Prevention

### Automatic Mechanism

From equation (1): `P(xᵢ) = Q(xᵢ)·u(xᵢ)`

If Q(xᵢ) → 0:
- Must have P(xᵢ) → 0 (since u(xᵢ) is finite)
- Enforced automatically by the equation!

**No need for explicit constraints.**

### Bounded Solution Values

Can add simple box constraints:
```
-M ≤ u(xi) ≤ M
```

This ensures solution stays bounded, preventing poles.

**Much simpler than constraining Q > δ everywhere!**

---

## Extensions

### 1. Higher-Order Equations

For 4th-order equation: `u⁽⁴⁾ = f(x)`

Add more differentiation:
```
P⁽³⁾ = Q⁽³⁾·u + 3Q''·u' + 3Q'·u'' + Q·u⁽³⁾
P⁽⁴⁾ = ... + Q·u⁽⁴⁾ = ... + Q·f
```

Still quadratic in unknowns!

### 2. Systems of ODEs

For `u₁'' = f₁(u₁, u₂)`, `u₂'' = f₂(u₁, u₂)`:

Define `P₁ = Q₁·u₁`, `P₂ = Q₂·u₂`

Get coupled quadratic system - still manageable!

### 3. Nonlinear ODEs

For `-u'' = f(x, u, u')`:

```
P'' = Q''·u + 2Q'·u' - Q·f(x, u, u')
```

If f is polynomial in u, u': still low-degree polynomial!

Example: `f(x, u, u') = u²` → cubic nonlinearity (not bad!)

---

## Implementation Recommendations

### Solver Choice

**For this quadratic formulation:**

1. **Levenberg-Marquardt** (first choice)
   - Designed for nonlinear least squares
   - Works well with quadratic structure
   - `scipy.optimize.least_squares(method='lm')`

2. **Trust-region methods**
   - Good for constrained problems
   - `scipy.optimize.minimize(method='trust-constr')`

3. **Sequential Quadratic Programming (SQP)**
   - Natural for quadratic constraints
   - `scipy.optimize.minimize(method='SLSQP')`

4. **Bilinear alternating**
   - Custom implementation
   - Alternately solve linear systems
   - Often fastest for well-conditioned problems

### Initial Guess Strategy

**Option 1: From polynomial solution**
```python
# Solve with finite differences
u_poly = solve_finite_difference(f, a, b, alpha, beta)

# Interpolate to get initial guess
for i, xi in enumerate(collocation_points):
    u_vals[i] = u_poly(xi)
    up_vals[i] = u_poly.deriv(1)(xi)

# Fit P and Q to match u_poly
P_init = fit_polynomial(u_poly)
Q_init = [1.0, 0.0, ...]  # Start with Q ≈ 1
```

**Option 2: From rational interpolation**
```python
# Use two-point Padé to interpolate boundary values
rational_init = TwoPointPade.from_endpoint_derivatives(...)
P_init = rational_init.numerator_coeffs()
Q_init = rational_init.denominator_coeffs()

# Evaluate at collocation points
for i, xi in enumerate(collocation_points):
    u_vals[i] = rational_init.evaluate(xi)
    up_vals[i] = rational_init.derivative(xi, 1)
```

---

## Comparison: All Three Formulations

### Residual Evaluation Complexity

**Quotient form:**
```python
uxx = (Pxx*Q**2 - 2*Px*Qx*Q - P*Qxx*Q + 2*P*Qx**2) / Q**3
```
- 1 division
- 5 multiplications
- 4 additions

**Cleared form:**
```python
lhs = Q**2*Pxx - 2*Q*Qx*Px + (2*Qx**2 - Q*Qxx)*P
rhs = -Q**3 * f
```
- 0 divisions
- 7 multiplications
- 4 additions

**Quadratic form:**
```python
res1 = P - Q*u
res2 = Px - Qx*u - Q*up
res3 = Pxx - Qxx*u - 2*Qx*up + Q*f
```
- 0 divisions
- 6 multiplications (across 3 equations)
- 6 additions

**Quadratic is simplest!**

### Degrees of Freedom

For [2/1] rational on [0,1] with k collocation points:

| Formulation | Unknowns | Equations | k required |
|-------------|----------|-----------|------------|
| Quotient | 4 (a₀,a₁,a₂,b₁) | 2 BC + 2k | any k ≥ 1 |
| Cleared | 4 | 2 BC + 2k | any k ≥ 1 |
| Quadratic | 4 + 2k | 2 BC + 3k | k = 2 (fixed) |

Quadratic has **more unknowns** but **lower degree**!

Trade-off: More variables but simpler structure.

---

## Numerical Experiments Needed

### Test Cases

1. **Smooth Poisson:** `-u'' = π²sin(πx)`
   - Compare convergence rates
   - Measure solve time

2. **Boundary layer:** `-εu'' + u' = 1` with ε = 0.01
   - Test stability near pole
   - Compare accuracy

3. **Nonlinear:** `-u'' = u²`
   - Assess advantage for nonlinear ODEs
   - Count iterations to convergence

### Metrics

- Solution accuracy: ||u_computed - u_exact||
- Convergence: Newton iterations vs mesh refinement
- Stability: Behavior as Q → 0
- Cost: Time per solve

---

## Conclusion

### The Key Insight

**By treating u and u' as explicit unknowns**, the system becomes:
- **Quadratic** instead of cubic
- **More structured** (bilinear)
- **Easier to solve** (specialized methods available)
- **Same stability** (no divisions)

### Recommendation

**For Gelfgren implementation:**

1. **Implement all three formulations** (quotient, cleared, quadratic)
2. **Default to quadratic** for general use
3. **Benchmark thoroughly** to confirm advantages
4. **Use bilinear alternating** for speed-critical applications

### Why This Matters

This quadratic formulation:
- ✅ Lowers polynomial degree (3 → 2)
- ✅ Enables specialized solvers (SQP, L-M)
- ✅ Provides clear physical interpretation (u values explicit)
- ✅ Natural for alternating optimization
- ✅ Easy to add regularization

**This could be the best way to implement rational collocation!**

---

**Document Status:** ✅ Complete analysis of quadratic formulation

**Next Steps:** Implement and benchmark against cleared form!
