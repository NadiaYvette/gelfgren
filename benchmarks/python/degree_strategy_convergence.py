#!/usr/bin/env python3
"""
Degree Progression Strategy Convergence Study

Compares three strategies:
1. BALANCED [n/m]: n ≈ 2m 
2. FIXED [n/2]: One complex conjugate pair
3. FIXED [n/4]: Two complex conjugate pairs

Tests on Smooth Poisson problem with systematic degree increases.
"""

import numpy as np
import matplotlib.pyplot as plt
import json
from dataclasses import dataclass, asdict
from typing import List, Tuple
import os

from rational_collocation_cleared import RationalCollocationCleared, QConstraintType


@dataclass
class StrategyResult:
    """Results for one degree configuration"""
    n: int
    m: int
    dof: int
    l2_error: float
    l_inf_error: float
    h1_seminorm_error: float
    iterations: int
    success: bool


@dataclass
class StrategyConvergence:
    """Complete convergence study for one strategy"""
    strategy_name: str
    description: str
    results: List[StrategyResult]


def smooth_poisson_forcing(x):
    """Smooth Poisson: -u'' = π² sin(πx)"""
    return np.pi ** 2 * np.sin(np.pi * x)


def smooth_poisson_exact(x):
    """Exact solution: u = sin(πx)"""
    return np.sin(np.pi * x)


def compute_errors(solver, result, a, b):
    """Compute L2, Linf, and H1 seminorm errors"""
    x_fine = np.linspace(a, b, 1000)
    u_approx = solver.evaluate_solution(result, x_fine)
    u_exact = smooth_poisson_exact(x_fine)
    
    h = x_fine[1] - x_fine[0]
    
    # L2 error
    l2_error = np.sqrt(h * np.sum((u_approx - u_exact) ** 2))
    
    # Linf error
    linf_error = np.max(np.abs(u_approx - u_exact))
    
    # H1 seminorm error
    du_approx = np.diff(u_approx) / h
    du_exact = np.diff(u_exact) / h
    h1_error = np.sqrt(h * np.sum((du_approx - du_exact) ** 2))
    
    return l2_error, linf_error, h1_error


def run_strategy_convergence(strategy_name, degree_pairs, max_dof=30):
    """Run convergence study for one strategy"""
    print(f"\n{strategy_name}")
    print("=" * 80)
    
    results = []
    for n, m in degree_pairs:
        dof = n + 1 + m
        if dof > max_dof:
            break
        
        try:
            solver = RationalCollocationCleared(
                f=smooth_poisson_forcing,
                a=0.0, b=1.0,
                alpha=0.0, beta=0.0,
                n=n, m=m,
                q_constraint=QConstraintType.ENDPOINT,
                constraint_epsilon=1e-4
            )
            
            result = solver.solve(method='trf', verbose=False)
            
            if result.success:
                l2, linf, h1 = compute_errors(solver, result, 0.0, 1.0)
                
                results.append(StrategyResult(
                    n=n, m=m, dof=dof,
                    l2_error=l2, l_inf_error=linf, h1_seminorm_error=h1,
                    iterations=result.iterations, success=True
                ))
                
                print(f"  [{n}/{m}] (DOF={dof:2d}): L2={l2:.6e}, L∞={linf:.6e}, H1={h1:.6e}, iter={result.iterations}")
            else:
                results.append(StrategyResult(
                    n=n, m=m, dof=dof,
                    l2_error=float('inf'), l_inf_error=float('inf'), h1_seminorm_error=float('inf'),
                    iterations=result.iterations, success=False
                ))
                print(f"  [{n}/{m}] (DOF={dof:2d}): FAILED - {result.message}")
        except Exception as e:
            print(f"  [{n}/{m}] (DOF={dof:2d}): ERROR - {str(e)[:50]}")
    
    return results


def main():
    print("Degree Progression Strategy Convergence Study")
    print("=" * 80)
    print("Problem: Smooth Poisson -u'' = π² sin(πx)")
    
    # Define strategies
    strategies = [
        ("Balanced [n/m]", "n ≈ 2m", 
         [(4, 2), (6, 3), (8, 4), (10, 5), (12, 6), (14, 7), (16, 8), (18, 9), (20, 10)]),
        
        ("Fixed [n/2]", "One complex conjugate pair",
         [(4, 2), (6, 2), (8, 2), (10, 2), (12, 2), (14, 2), (16, 2), (18, 2), (20, 2), (22, 2), (24, 2)]),
        
        ("Fixed [n/4]", "Two complex conjugate pairs",
         [(8, 4), (10, 4), (12, 4), (14, 4), (16, 4), (18, 4), (20, 4), (22, 4), (24, 4), (26, 4)])
    ]
    
    # Run convergence studies
    all_results = []
    for name, desc, degrees in strategies:
        results = run_strategy_convergence(name, degrees, max_dof=35)
        all_results.append(StrategyConvergence(name, desc, results))
    
    # Save results
    output_dir = os.path.join(os.path.dirname(__file__), '../data')
    os.makedirs(output_dir, exist_ok=True)
    
    output_file = os.path.join(output_dir, 'degree_strategy_convergence.json')
    
    # Convert to serializable format
    data = {
        'problem': 'Smooth Poisson',
        'strategies': [
            {
                'name': strat.strategy_name,
                'description': strat.description,
                'results': [asdict(r) for r in strat.results]
            }
            for strat in all_results
        ]
    }
    
    with open(output_file, 'w') as f:
        json.dump(data, f, indent=2)
    
    print(f"\n{'=' * 80}")
    print(f"Results saved to {output_file}")
    
    # Generate plots
    plot_strategy_comparison(all_results)
    
    print("=" * 80)
    print("Convergence study complete!")


def plot_strategy_comparison(all_results):
    """Generate comparison plots"""
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 10))
    
    colors = ['blue', 'red', 'green']
    markers = ['o', 's', '^']
    
    for idx, strat in enumerate(all_results):
        # Extract successful results
        dofs = [r.dof for r in strat.results if r.success]
        l2_errors = [r.l2_error for r in strat.results if r.success]
        linf_errors = [r.l_inf_error for r in strat.results if r.success]
        
        if not dofs:
            continue
        
        # Plot 1: L2 error vs DOF
        ax1.loglog(dofs, l2_errors, f'{markers[idx]}-', 
                   label=strat.strategy_name, color=colors[idx], linewidth=2, markersize=8)
        
        # Plot 2: Linf error vs DOF
        ax2.loglog(dofs, linf_errors, f'{markers[idx]}-',
                   label=strat.strategy_name, color=colors[idx], linewidth=2, markersize=8)
        
        # Plot 3: L2 error vs configuration index
        ax3.semilogy(range(len(l2_errors)), l2_errors, f'{markers[idx]}-',
                     label=strat.strategy_name, color=colors[idx], linewidth=2, markersize=8)
        
        # Plot 4: Convergence rate
        if len(l2_errors) > 1:
            rates = []
            dof_mids = []
            for i in range(1, len(l2_errors)):
                if l2_errors[i] > 0 and l2_errors[i-1] > 0:
                    rate = np.log(l2_errors[i-1] / l2_errors[i]) / np.log(dofs[i] / dofs[i-1])
                    rates.append(rate)
                    dof_mids.append((dofs[i] + dofs[i-1]) / 2)
            
            if rates:
                ax4.plot(dof_mids, rates, f'{markers[idx]}-',
                         label=strat.strategy_name, color=colors[idx], linewidth=2, markersize=8)
    
    # Configure axes
    ax1.set_xlabel('Degrees of Freedom', fontsize=12)
    ax1.set_ylabel('L² Error', fontsize=12)
    ax1.set_title('L² Error vs DOF', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    
    ax2.set_xlabel('Degrees of Freedom', fontsize=12)
    ax2.set_ylabel('L∞ Error', fontsize=12)
    ax2.set_title('L∞ Error vs DOF', fontsize=14, fontweight='bold')
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)
    
    ax3.set_xlabel('Configuration Index', fontsize=12)
    ax3.set_ylabel('L² Error (log scale)', fontsize=12)
    ax3.set_title('Error Progression', fontsize=14, fontweight='bold')
    ax3.legend(fontsize=10)
    ax3.grid(True, alpha=0.3)
    
    ax4.set_xlabel('Degrees of Freedom', fontsize=12)
    ax4.set_ylabel('Convergence Rate α', fontsize=12)
    ax4.set_title('Convergence Rate (error ~ DOF^-α)', fontsize=14, fontweight='bold')
    ax4.legend(fontsize=10)
    ax4.grid(True, alpha=0.3)
    ax4.axhline(y=0, color='k', linestyle='--', alpha=0.3)
    
    plt.tight_layout()
    
    # Save figure
    figures_dir = os.path.join(os.path.dirname(__file__), '../reports/figures')
    os.makedirs(figures_dir, exist_ok=True)
    output_file = os.path.join(figures_dir, 'degree_strategy_comparison.pdf')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Figure saved to {output_file}")
    plt.close()


if __name__ == '__main__':
    main()
