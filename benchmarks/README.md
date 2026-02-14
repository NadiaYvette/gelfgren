# Gelfgren Benchmarks

> **✅ Python bindings are now available!** The special functions benchmark (`special_function_convergence.py`) uses the real Gelfgren implementation and demonstrates significant improvements for rational approximants over polynomial splines.
>
> **Note:** The BVP benchmark (`bvp_convergence.py`) still uses placeholder code as it requires implementing a full rational-based BVP solver, which is beyond the scope of an interpolation library.

Comprehensive benchmark suite comparing piecewise rational approximants against polynomial splines for:

1. **Boundary Value Problems (BVPs)**: 1D Poisson equation variants
2. **Special Function Approximations**: Exponential, trigonometric, error function, Bessel functions, etc.

## Overview

These benchmarks test the hypothesis that **piecewise rational approximants can achieve comparable accuracy to polynomial splines with coarser meshes**, thanks to their greater flexibility in representing functions with poles, rapid variation, and other challenging features.

## Structure

```
benchmarks/
├── Cargo.toml                      # Rust benchmark package
├── lib.rs                          # Main benchmark library
├── bvp/                            # Boundary value problems
│   └── poisson_1d.rs              # 1D Poisson equation test cases
├── special_functions/              # Special function approximations
│   └── approximations.rs          # Various special functions
├── python/                         # Python convergence studies
│   ├── bvp_convergence.py         # BVP convergence analysis
│   ├── special_function_convergence.py  # Special function analysis
│   └── generate_latex_report.py   # LaTeX report generator
├── data/                           # Generated benchmark data (JSON)
│   ├── *.json                     # BVP results
│   └── special_functions/         # Special function results
│       └── *.json
└── reports/                        # Generated reports
    ├── figures/                   # Convergence plots (PDF)
    │   ├── *.pdf                  # BVP plots
    │   └── special_functions/     # Special function plots
    │       └── *.pdf
    └── latex/                     # LaTeX source and PDFs
        ├── comprehensive_benchmark_report.tex
        └── comprehensive_benchmark_report.pdf
```

## Running the Benchmarks

### Prerequisites

- **Rust**: 1.70+ with Cargo
- **Python**: 3.9+ with NumPy, SciPy, Matplotlib
- **LaTeX** (optional): For PDF report generation (pdflatex)

### Python Environment Setup

```bash
# Install Python dependencies
pip install numpy scipy matplotlib

# Build Gelfgren Python bindings
cd ../bindings/python
pip install maturin
maturin develop
cd ../../benchmarks
```

### Run All Benchmarks

```bash
# Run BVP convergence studies
python python/bvp_convergence.py

# Run special function approximation studies
python python/special_function_convergence.py

# Generate comprehensive LaTeX report
python python/generate_latex_report.py --mode comprehensive

# Compile to PDF (requires pdflatex)
cd reports/latex
pdflatex comprehensive_benchmark_report.tex
pdflatex comprehensive_benchmark_report.tex  # Second pass for references
```

### Individual Components

**BVP benchmarks only:**
```bash
python python/bvp_convergence.py
python python/generate_latex_report.py --mode bvp
```

**Special functions only:**
```bash
python python/special_function_convergence.py
```

## Benchmark Problems

### Boundary Value Problems

All BVPs solve: $-u''(x) = f(x)$ on $[0, 1]$ with $u(0) = u(1) = 0$

1. **Smooth Poisson**: $f(x) = \pi^2 \sin(\pi x)$
   - Exact solution: $u(x) = \sin(\pi x)$
   - Tests smooth, well-behaved case

2. **Discontinuous Poisson**: Piecewise constant forcing
   - Tests handling of discontinuities
   - Non-smooth solution

3. **Oscillatory Poisson**: $f(x) = (10\pi)^2 \sin(10\pi x)$
   - High-frequency oscillations
   - Challenges for coarse meshes

### Special Functions

1. **Exponential**: $e^x$ on $[-1, 1]$ and $[-5, 5]$
   - Rapid growth
   - Well-suited to rational approximation

2. **Trigonometric**: $\sin(x)$, $\cos(x)$, $\tan(x)$
   - Periodic and smooth
   - Tangent has poles

3. **Error Function**: $\text{erf}(x)$ on $[-3, 3]$
   - No closed form
   - Smooth sigmoid shape

4. **Bessel Function**: $J_0(x)$ on $[0, 10]$
   - Oscillatory decay
   - Arises in wave propagation

5. **Logarithm**: $\log(1+x)$ on $[0, e-1]$
   - Pole at $x = -1$
   - Slow variation

6. **Runge's Function**: $\frac{1}{1 + 25x^2}$ on $[-1, 1]$
   - Classic example of Runge's phenomenon
   - Polynomial interpolation fails badly
   - Rational approximants should excel

## Error Metrics

For each test case and mesh size, we compute:

- **L² error**: $\|u - u_h\|_{L^2} = \sqrt{\int_a^b |u(x) - u_h(x)|^2 \, dx}$
- **L∞ error**: $\|u - u_h\|_{L^\infty} = \max_{x \in [a,b]} |u(x) - u_h(x)|$
- **H¹ seminorm error**: $|u - u_h|_{H^1} = \|u' - u_h'\|_{L^2}$
- **Relative errors**: Normalized by exact solution norm

## Convergence Analysis

For each method (polynomial splines vs piecewise rational approximants), we:

1. **Vary mesh size**: $n = 4, 8, 16, 32, 64, 128$ intervals
2. **Compute error metrics** at each refinement level
3. **Calculate convergence rates**: $\alpha$ where error $\sim h^\alpha$
4. **Compare DOF efficiency**: Error vs degrees of freedom

### Degrees of Freedom

- **Polynomial (cubic spline)**: DOF = $n + 3$ (for $n$ intervals)
- **Rational ([2/1] Padé per interval)**: DOF = $4n$ (accounting for normalization constraint)
- **Rational ([3/2] Padé per interval)**: DOF = $6n$ (same DOF per interval as quintic polynomials)

**Note on [2/2]:** The [2/2] rational is not achievable with two-point Padé! The constraint $n+m+1=2p$ (must be even) means [2/2] with $2+2+1=5$ is impossible. The closest balanced rationals are [2/1] or [1/2] with 4 DOF per interval.

Rational approximants use more DOF per interval but may achieve target accuracy with fewer intervals (coarser mesh). The [3/2] rational has the same DOF per interval as quintic polynomials but with potentially better approximation properties for certain function classes.

## Key Questions

The benchmarks address:

1. **Do rational approximants achieve comparable accuracy with coarser meshes?**
   - Compare error for same mesh size
   - Find equivalent meshes for target accuracy

2. **What is the convergence rate?**
   - Polynomial splines: Expected $O(h^4)$ for smooth problems
   - Rational approximants: May achieve higher rates

3. **Which problems benefit most from rational approximants?**
   - Functions with poles or near-poles
   - Rapidly varying functions
   - Problems where polynomial interpolation struggles

## Output

### Data Files (JSON)

Located in `data/`:
- Convergence data for each problem
- Error metrics at each mesh size
- Computed convergence rates

### Figures (PDF)

Located in `reports/figures/`:
- Error vs mesh size (log-log)
- Error vs DOF (efficiency comparison)
- Convergence rate plots

### LaTeX Report (PDF)

Located in `reports/latex/comprehensive_benchmark_report.pdf`:
- Mathematical problem formulations with equations
- Convergence tables
- Embedded figures
- Theoretical analysis
- Conclusions and recommendations

## Theoretical Background

### Polynomial Splines

For smooth functions $u \in C^{n+1}[a,b]$, cubic spline interpolation achieves:

$$\|u - s\|_{L^2} \leq C h^4 \|u^{(4)}\|_{L^2}$$

where $h$ is the mesh size and $s$ is the interpolating spline.

### Rational Approximants

Padé approximants $[m/n]$ can represent functions with poles exactly (if the pole structure matches). For meromorphic functions, they may achieve exponential convergence.

### Runge's Phenomenon

High-degree polynomial interpolation at equidistant points can exhibit wild oscillations near boundaries. Rational approximants avoid this issue.

## Reproducibility

All code is deterministic. Benchmark results should be reproducible given:
- Same Rust/Python versions
- Same random number generator seeds (none used currently)
- Same numerical libraries

## Citation

If you use these benchmarks in your research, please cite:

```
@software{gelfgren_benchmarks,
  title = {Gelfgren Benchmarks: Piecewise Rational vs Polynomial Approximation},
  author = {Chambers, Nadia and Claude AI},
  year = {2025},
  url = {https://github.com/yourusername/gelfgren}
}
```

## License

Licensed under either of:
- Apache License, Version 2.0 ([LICENSE-APACHE](../LICENSE-APACHE))
- MIT License ([LICENSE-MIT](../LICENSE-MIT))

at your option.

## Contributing

Contributions welcome! Areas for expansion:
- Additional BVP types (2D problems, eigenvalue problems)
- More special functions (hypergeometric, elliptic integrals)
- Adaptive mesh refinement studies
- Performance benchmarks (timing, memory)

## References

1. Gelfgren, J. (1975). "Piecewise Rational Interpolation". *BIT Numerical Mathematics* 15:382–393.
2. de Boor, C. (2001). *A Practical Guide to Splines*. Springer.
3. Baker, G.A. and Graves-Morris, P. (1996). *Padé Approximants*. Cambridge University Press.
4. Traub, J.F. (1964). "On Lagrange-Hermite Interpolation". *SIAM Journal on Numerical Analysis*.
