# Extended Benchmark Analysis Summary

## Best Method for Each Problem

| Problem | Best Method | Error | Convergence Rate |
|---------|-------------|-------|------------------|
| AdvDiff_Mild_eps0.1 | Polynomial FD 4th | 2.70e-03 | 1.77 |
| AdvDiff_Sharp_eps0.001 | Polynomial FD 4th | 7.77e-02 | 2.25 |
| Helmholtz_Exp_k4 | Rational Cleared | 2.03e-07 | 30.24 |
| Helmholtz_Sin_k1 | Rational Cleared | 1.84e-12 | 21.25 |
| ReacDiff_Exp_c10 | Rational Cleared | 4.52e-13 | 19.75 |
| ReacDiff_Polynomial_c4 | Rational Cleared | 4.23e-17 | -9.19 |
| VarCoeff_Discontinuous | Rational Cleared | 3.90e-01 | 0.06 |
| VarCoeff_Poly | Polynomial FD 2nd | 3.52e-01 | 0.00 |

## Convergence Rate Summary

| Problem | 2nd FD | 4th FD | 6th FD | Cheby | Leg | Rational |
|---------|--------|--------|--------|-------|-----|----------|
| AdvDiff_Mild_eps0.1 | 0.00 | 1.77 | N/A | -0.00 | -0.00 | 4.43 |
| AdvDiff_Sharp_eps0.001 | 0.00 | 2.25 | N/A | 0.76 | 1.12 | 0.32 |
| Helmholtz_Exp_k4 | -0.01 | 1.88 | N/A | 2.08 | 2.04 | 30.24 |
| Helmholtz_Sin_k1 | 0.00 | 1.79 | N/A | 2.08 | 2.04 | 21.25 |
| ReacDiff_Exp_c10 | -0.00 | 1.97 | N/A | 2.08 | 2.04 | 19.75 |
| ReacDiff_Polynomial_c4 | 0.01 | 2.34 | N/A | 2.08 | 2.04 | -9.19 |
| VarCoeff_Discontinuous | -0.00 | N/A | N/A | N/A | N/A | 0.06 |
| VarCoeff_Poly | 0.00 | N/A | N/A | N/A | N/A | -0.15 |
