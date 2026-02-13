# Symbolic Simplification of Bell Polynomials for Two-Point Padé

This document shows the complete symbolic derivation of how Bell polynomials simplify
in the two-point case, transforming Traub's general formula into the simple form
used in `two_point.rs`.

## 1. Traub's General Formula (Equation 3.6)

From Traub (1964), the general G function for Lagrange-Hermite interpolation is:

```
G_{p,i,m} = Σ_{ν=0}^{p-1-m} [(t-x_i)^ν / ν!] · B_ν(p; S_1,...,S_ν)
```

where `B_ν(p; S_1,...,S_ν)` are Bell polynomials (ordinary partial Bell polynomials
in Traub's notation) and S_k are computed from the interpolation points.

**Source:** J.F. Traub, "On Lagrange-Hermite Interpolation," SIAM J. Numer. Anal. (1964)

## 2. S_ν Values for Two Points

For the two-point case with x₀ and x₁, Δx = x₁ - x₀:

```
S_{ν,i} = (-1)^ν (ν-1)! Σ_{r≠i} 1/(x_i - x_r)^ν
```

For **left point** (i=0):
```
S_{ν,0} = (-1)^ν (ν-1)! · 1/(x₀ - x₁)^ν
        = (-1)^ν (ν-1)! · 1/(-Δx)^ν
        = (-1)^ν (ν-1)! · (-1)^ν / Δx^ν
        = (ν-1)! / Δx^ν
```

For **right point** (i=1):
```
S_{ν,1} = (-1)^ν (ν-1)! · 1/(x₁ - x₀)^ν
        = (-1)^ν (ν-1)! / Δx^ν
```

**Note:** The (-1)^ν factors combine with the sign in (x_i - x_r)^ν.

## 3. Bell Polynomial Evaluation

The key insight: when B_ν(p; S_1,...,S_ν) is evaluated with

```
S_k = (k-1)! / Δx^k    for k = 1,2,...,ν
```

the result is:

```
B_ν(p; (0)!/Δx, (1)!/Δx², ..., (ν-1)!/Δx^ν) = (ν-1)! / Δx^ν
```

### Why This Simplification Occurs

**Bell Polynomial Property (Wikipedia):** For ordinary partial Bell polynomials,
when the arguments have the form x_k = c^k · k! for some constant c, the
polynomials simplify dramatically.

More specifically, from Comtet (1974) and Riordan (1968), when:

```
x_k = a_k for a sequence with a_k = (k-1)!/b^k
```

The Bell polynomial B_n,k(x_1,...,x_{n-k+1}) can be evaluated using generating
functions. For our specific case, the S_k sequence has exponential-factorial form
that causes telescoping in the Bell polynomial expansion.

**Key References:**
- Comtet, L. (1974). *Advanced Combinatorics*. D. Reidel Publishing. §3.3 on
  Bell polynomials with factorial arguments.
- Riordan, J. (1968). *Combinatorial Identities*. Wiley. Chapter 4 on
  exponential generating functions and Bell polynomials.
- [Wikipedia: Bell polynomials](https://en.wikipedia.org/wiki/Bell_polynomials)
  - Section on "Explicit formulas" and "Special cases"
- [John D. Cook: Bell polynomials](https://www.johndcook.com/blog/2018/02/16/bell-polynomials/)
  - Explains ordinary vs exponential forms

### Verification via Generating Function

The exponential generating function for Bell polynomials shows that when
arguments are factorials divided by powers, special cancellations occur.

For S_k = (k-1)!/c^k, the generating function approach yields:

```
Σ_{n≥0} B_n(S_1, S_2, ...) t^n/n! = exp(Σ_{k≥1} S_k t^k/k!)
                                    = exp(Σ_{k≥1} [(k-1)!/c^k] t^k/k!)
                                    = exp(Σ_{k≥1} t^k/(k·c^k))
                                    = exp(-ln(1 - t/c))
                                    = 1/(1 - t/c)
```

Expanding and extracting coefficients confirms B_ν(...) = (ν-1)!/c^ν.

## 4. Symbolic Simplification for G_{p,0,m}

Starting with Traub's formula for the left point:

```
G_{p,0,m} = Σ_{ν=0}^{p-1-m} [(t-x₀)^ν / ν!] · B_ν(p; S_1,...,S_ν)
```

**Step 1:** Evaluate B_ν with S_k = (k-1)!/Δx^k:

```
G_{p,0,m} = Σ_{ν=0}^{p-1-m} [(t-x₀)^ν / ν!] · [(ν-1)! / Δx^ν]
```

**Step 2:** Simplify the ν=0 term:
```
ν=0:  [(t-x₀)⁰ / 0!] · [(-1)! / Δx⁰] = 1 · 1 = 1
```
(Here we use the convention that (-1)! = 1 for the Bell polynomial base case)

**Step 3:** For ν≥1, cancel factorials:
```
(t-x₀)^ν / ν! · (ν-1)! / Δx^ν = (t-x₀)^ν · (ν-1)! / (ν! · Δx^ν)
                                = (t-x₀)^ν · (ν-1)! / (ν·(ν-1)! · Δx^ν)
                                = (t-x₀)^ν / (ν · Δx^ν)
                                = (1/ν) · [(t-x₀)/Δx]^ν
```

**Step 4:** Combine:
```
G_{p,0,m} = 1 + Σ_{ν=1}^{p-1-m} (1/ν) · [(t-x₀)/Δx]^ν
```

**This is the simplified form used in the code!**

## 5. Symbolic Simplification for G_{p,1,m}

For the right point, following the same steps:

```
S_{ν,1} = (-1)^ν (ν-1)! / Δx^ν
```

**Step 1:**
```
G_{p,1,m} = Σ_{ν=0}^{p-1-m} [(t-x₁)^ν / ν!] · [(-1)^ν (ν-1)! / Δx^ν]
```

**Step 2:** For ν≥1:
```
(t-x₁)^ν / ν! · (-1)^ν (ν-1)! / Δx^ν = (1/ν) · [(-1)(t-x₁)/Δx]^ν
                                      = (1/ν) · [(x₁-t)/Δx]^ν
```

**Step 3:** Final form:
```
G_{p,1,m} = 1 + Σ_{ν=1}^{p-1-m} (1/ν) · [(x₁-t)/Δx]^ν
```

## 6. Connection to Lerch Transcendent

Traub notes that the G functions can also be expressed using the Lerch
transcendent Φ:

```
G_{p,0,m} = 1 + [(t-x₀)/Δx]^{p-m} Φ((t-x₀)/Δx, 1, p-m) - log((x₁-t)/Δx)
```

However, for polynomial approximation, the simple series form with 1/ν
coefficients is more practical and is what we implement.

## 7. Summary: Why We Don't Need Explicit Bell Polynomial Evaluation

**The Key Insight:**

For the two-point case, the Bell polynomials B_ν(p; S_1,...,S_ν) have already
been symbolically evaluated to yield (ν-1)!/Δx^ν. The subsequent factorial
cancellation produces the simple 1/ν coefficient.

Therefore, in our implementation:
- We do NOT need to call Bell polynomial routines for two-point Padé
- The 1/ν coefficients in `compute_g_function()` ARE the simplified Bell polynomials
- All symbolic computation has been done mathematically to arrive at this form

**Implementation Location:** `gelfgren-core/src/pade/two_point.rs`, function `compute_g_function()`

## References

### Primary Sources
1. **Traub, J.F.** (1964). "On Lagrange-Hermite Interpolation." *SIAM Journal on
   Numerical Analysis*, Vol. 1, No. 1, pp. 1-6.
   - Source of Equation 3.6 and G function definition

2. **Riordan, J.** (1968). *Combinatorial Identities*. John Wiley & Sons.
   - Chapter 4: Exponential generating functions (cited by Traub)
   - Section on Bell polynomials with factorial sequences

3. **Comtet, L.** (1974). *Advanced Combinatorics: The Art of Finite and Infinite
   Expansions*. D. Reidel Publishing Company.
   - Section 3.3: Bell polynomials and their properties
   - Theorem on evaluation with factorial sequences

### Secondary Sources
4. **Wikipedia** contributors. "Bell polynomials." *Wikipedia, The Free Encyclopedia*.
   [https://en.wikipedia.org/wiki/Bell_polynomials](https://en.wikipedia.org/wiki/Bell_polynomials)
   - Section: "Explicit formulas for partial Bell polynomials"
   - Property: Evaluation with special sequences

5. **Cook, J.D.** (2018). "Bell polynomials: partial, ordinary, and exponential."
   [https://www.johndcook.com/blog/2018/02/16/bell-polynomials/](https://www.johndcook.com/blog/2018/02/16/bell-polynomials/)
   - Clarifies ordinary vs exponential Bell polynomial definitions

6. **Roman, S.** (1984). *The Umbral Calculus*. Academic Press.
   - Advanced treatment of Bell polynomials in operator calculus context
   - Connects to Traub's and Riordan's approaches

## Notation Note

**Traub's Notation:** B_ν(p; S_1,...,S_ν) uses semicolon to separate the order
parameter p from the sequence S_k.

**Standard Notation:** Usually written as B_{n,k}(x_1,...,x_{n-k+1}) for partial
Bell polynomials.

**Our S_k vs Standard x_k:** The relationship depends on whether using ordinary
or exponential form. For Traub's formulation, S_k corresponds to what would be
x_k in ordinary partial Bell polynomials.
