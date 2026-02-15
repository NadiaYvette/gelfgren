# Special Function Approximation Benchmark Results

## Summary

Comprehensive benchmarks comparing **piecewise polynomial splines** vs **piecewise rational approximants** for approximating 16 special functions across 6 mesh refinement levels (4, 8, 16, 32, 64, 128 intervals).

## Key Findings

### 1. Polynomial Collocation Dominates (10/16 wins)

**Winner functions:**
- Sine, Cosine (both ~15× better at 128 intervals)
- Bessel J₁, J₅, J₁₀ (12-29× better)
- Mathieu ce₀, se₂ (100-200× better, rationals stagnate)
- Jacobi sn, Lemniscate sl (17× better)
- Airy Ai oscillatory (17× better)
- **Runge's function** (27× better - surprising!)

**Convergence:** Consistent α ≈ 4.4 across all functions

### 2. Rational Methods: Limited Advantage (1/16 wins)

**Winner functions:**
- Logarithm log(1+x) only (33× better)

**Comparable (5/16):**
- Exponential e^x, Error function, Cosine, Jacobi cn, Airy Bi

**Critical failures:**
- Mathieu ce₀: Convergence stagnates (α = 1.4 vs polynomial α = 4.4)
- High-order Bessel: Struggle with varying frequency
- Runge's function: Piecewise polynomials avoid global pathology

### 3. Surprising Results

#### Runge's Function Contradiction
- **Expected:** Rationals excel (function naturally rational: 1/(1+25x²))
- **Actual:** Polynomials win decisively (4.5×10⁻⁷ vs 1.2×10⁻⁵)
- **Reason:** Piecewise approach avoids global Runge phenomenon

#### Mathieu Convergence Stagnation
- Polynomial: Achieves machine precision (10⁻¹⁴)
- Rational: Stalls at 10⁻¹² with α = 1.4
- **Interpretation:** Complex oscillatory structure challenges rational basis

#### High-Order Bessel Failures
- J₁₀: Polynomial 29× better despite smooth oscillations
- Varying oscillatory frequency defeats rational basis uniformity

### 4. Efficiency Comparison

At 128 intervals (finest mesh):
- **Polynomial DOF:** 131 (2n + 3)
- **Rational DOF:** 512 (4n + 4)
- **Advantage:** Polynomials use 4× fewer degrees of freedom

### 5. Contrast with BVP Solving

| Problem Type | Polynomial | Rational | Why? |
|--------------|-----------|----------|------|
| **Smooth BVPs** (Chapters 11-15) | Good | **Excellent** (machine precision) | Cleared form provides spectral discretization |
| **Function Approximation** (Chapter 5) | **Excellent** | Mediocre | Direct interpolation has no such advantage |

**Key Insight:** Rational collocation's advantage is in the **cleared form formulation for differential operators**, not in direct function representation.

## Recommendations

### For Function Approximation

**Use Polynomial Splines (cubic Hermite or higher) as default:**
- Robust: Consistent convergence (α ≈ 4.4)
- Efficient: 4× fewer DOF
- Predictable: No stagnation or divergence
- Accurate: Achieves 10⁻⁷ to 10⁻¹⁴ on all smooth functions

**Consider Rationals only for:**
- Logarithmic functions (demonstrated advantage)
- Functions with known exponential character where Padé theory applies
- Specific literature recommendations

**Avoid Rationals for:**
- Oscillatory functions (Bessel, Mathieu, Airy oscillatory)
- High-order functions (convergence rate unpredictable)
- General-purpose libraries (polynomials more robust)

### For BVP Solving

**Use Rational Collocation when:**
- Smooth problems where machine precision critical
- Spectral convergence needed
- Natural rational structure (exponential, reaction-diffusion)

**Use Polynomial/FD/Spectral when:**
- Non-smooth problems (boundary layers, discontinuities)
- Robustness priority over extreme accuracy
- Oscillatory forcing or solutions

## Technical Details

- **Polynomial method:** Piecewise cubic Hermite splines, C¹ continuous
- **Rational method:** P(x)/Q(x) with cubic numerator, linear denominator
- **Error metrics:** L² norm, L^∞ norm, H¹ seminorm, relative errors
- **Convergence rates:** Log-log regression, α where error ~ h^α
- **Mesh sequence:** [4, 8, 16, 32, 64, 128] intervals

## Files Generated

### Data Files (16 functions)
- `benchmarks/data/special_functions/*.json` - Convergence data

### Figures (16 functions)
- `benchmarks/reports/figures/special_functions/*.pdf` - Convergence plots

### Report
- Chapter 5 of `comprehensive_benchmark_report.tex` - Complete analysis
- Updated to 121 pages with actual results

## Broader Implications

1. **Piecewise > Global:** Piecewise polynomial approach avoids classical interpolation pathologies (Runge phenomenon)

2. **Method-Problem Match:** Rational methods are **specialist tools** requiring:
   - Guaranteed smoothness
   - Need for extreme accuracy or natural rational representation
   - Acceptance of higher computational cost

3. **Unified Smoothness Criterion:** Both BVP and function approximation show same pattern:
   - Smooth exponentials → Rationals competitive
   - Oscillatory/varying → Polynomials superior
   - Sharp gradients → Rationals catastrophic (BVPs only)

4. **Different Mechanisms:**
   - BVP solving: Rational advantage from cleared form operator discretization
   - Function approximation: No such advantage; polynomials inherently strong

## Reproducibility

Run benchmarks:
```bash
cd benchmarks/python
python special_function_convergence.py
```

Generate analysis:
```bash
python analyze_special_functions.py
```

Regenerate report:
```bash
cd benchmarks/reports/latex
pdflatex comprehensive_benchmark_report.tex
```
