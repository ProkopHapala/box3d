#!/usr/bin/env python3
"""
Didactic visualization of the hard/soft potential splitting used in RigidAtomFF.

The interatomic potential is:
    E(r) = q_i*q_j / r + C / r^n

where:
  - q_i*q_j / r  is the Coulomb term (attractive when q_i*q_j < 0)
  - C / r^n      is the short-range repulsion (C > 0, n = 4, 6, or 8)

Two fitting strategies for the hard parabolic spring:

  1. "r0" — Fit at the equilibrium r0 (C2-continuous, but parabola is too soft away from r0)
  2. "rfit" — Fit at r_fit = alpha * r0 (e.g. 0.9*r0), where the potential is already repulsive.
     The curvature there is higher, so the parabola better approximates the repulsive wall.
     The hard spring only activates for r < r_fit (shorter cutoff), and the soft potential
     handles the attractive well + transition region.

At r_fit, the parabola matches:
  - Value:  E_hard(r_fit) = E(r_fit)
  - Force:  F_hard(r_fit) = F(r_fit)  (nonzero — the spring has a bias force)
  - Curvature: k_hard = d²E/dr² at r_fit

The parabola with bias is:
  E_hard(r) = E(r_fit) + E'(r_fit)*(r - r_fit) + 0.5 * k * (r - r_fit)^2   for r < r_fit
  E_hard(r) = 0                                                            for r >= r_fit

In the SI solver, the linear term E'(r_fit) becomes a constant bias velocity,
and k becomes the spring stiffness with rest length r_fit.

Usage:
    python3 plot_hard_soft_split.py
    python3 plot_hard_soft_split.py --qprod -1.0 --C 1.0 --n 6 --rfit 0.9
    python3 plot_hard_soft_split.py --noPlot --saveFig hard_soft_split.png
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse


# ======================== Physics ========================

def E_full(r, qprod, C, n):
    """Full interatomic potential: Coulomb + power-law repulsion."""
    return qprod / r + C / r**n


def dE_full(r, qprod, C, n):
    """Derivative dE/dr of the full potential."""
    return -qprod / r**2 - n * C / r**(n + 1)


def d2E_full(r, qprod, C, n):
    """Second derivative d²E/dr² of the full potential."""
    return 2 * qprod / r**3 + n * (n + 1) * C / r**(n + 2)


def F_full(r, qprod, C, n):
    """Force F = -dE/dr (positive = repulsive, negative = attractive)."""
    return -dE_full(r, qprod, C, n)


def compute_r0(qprod, C, n):
    """
    Equilibrium distance where dE/dr = 0.
    Requires qprod < 0 (attractive Coulomb) for a minimum to exist.
    r0 = (n * C / |qprod|)^(1/(n-1))
    """
    assert qprod < 0, "Need qprod < 0 (attractive Coulomb) for a minimum to exist"
    return (n * C / abs(qprod))**(1.0 / (n - 1))


# ======================== Hard potential: fit at r0 (original) ========================

def E_hard_r0(r, qprod, C, n, r0):
    """
    Hard potential fit at r0: simple unilateral spring, zero force at r0.
    E_hard(r) = 0.5 * k * (r0 - r)^2  if r < r0, else 0.
    k = d²E/dr² at r0 = n*(n-1)*C / r0^(n+2)
    """
    k = n * (n - 1) * C / r0**(n + 2)
    out = np.zeros_like(r)
    mask = r < r0
    out[mask] = 0.5 * k * (r0 - r[mask])**2
    return out


def F_hard_r0(r, qprod, C, n, r0):
    """Hard force for r0 fit: F = k*(r0 - r) for r < r0, else 0."""
    k = n * (n - 1) * C / r0**(n + 2)
    out = np.zeros_like(r)
    mask = r < r0
    out[mask] = k * (r0 - r[mask])
    return out


# ======================== Hard potential: fit at r_fit (repulsive region) ========================

def E_hard_rfit(r, qprod, C, n, r_fit):
    """
    Hard potential fit at r_fit: parabola matching value, force, and curvature at r_fit.
    E_hard(r) = E(r_fit) + E'(r_fit)*(r - r_fit) + 0.5 * k * (r - r_fit)^2  if r < r_fit, else 0.
    k = d²E/dr² at r_fit.

    The linear term (E'(r_fit) != 0) provides a constant bias force so the parabola
    is tangent to E at r_fit, not just matching curvature.
    """
    r_fit_arr = np.array([r_fit])
    E_fit = E_full(r_fit_arr, qprod, C, n)[0]
    dE_fit = dE_full(r_fit_arr, qprod, C, n)[0]  # dE/dr at r_fit (negative => force is positive)
    k = d2E_full(r_fit_arr, qprod, C, n)[0]

    out = np.zeros_like(r)
    mask = r < r_fit
    dr = r[mask] - r_fit
    out[mask] = E_fit + dE_fit * dr + 0.5 * k * dr**2
    return out


def F_hard_rfit(r, qprod, C, n, r_fit):
    """
    Hard force for r_fit: F = -(dE_hard/dr) = -(E'(r_fit) + k*(r - r_fit))  for r < r_fit, else 0.
    At r_fit, F_hard = -E'(r_fit) = F_full(r_fit) — matches the real force.
    """
    r_fit_arr = np.array([r_fit])
    dE_fit = dE_full(r_fit_arr, qprod, C, n)[0]
    k = d2E_full(r_fit_arr, qprod, C, n)[0]

    out = np.zeros_like(r)
    mask = r < r_fit
    out[mask] = -(dE_fit + k * (r[mask] - r_fit))
    return out


# ======================== Soft residuals ========================

def E_soft_r0(r, qprod, C, n, r0):
    """Soft residual for r0 fit: E_soft = E - E_hard_r0."""
    return E_full(r, qprod, C, n) - E_hard_r0(r, qprod, C, n, r0)


def F_soft_r0(r, qprod, C, n, r0):
    """Soft force for r0 fit: F_soft = F - F_hard_r0."""
    return F_full(r, qprod, C, n) - F_hard_r0(r, qprod, C, n, r0)


def E_soft_rfit(r, qprod, C, n, r_fit):
    """Soft residual for r_fit: E_soft = E - E_hard_rfit."""
    return E_full(r, qprod, C, n) - E_hard_rfit(r, qprod, C, n, r_fit)


def F_soft_rfit(r, qprod, C, n, r_fit):
    """Soft force for r_fit: F_soft = F - F_hard_rfit."""
    return F_full(r, qprod, C, n) - F_hard_rfit(r, qprod, C, n, r_fit)


# ======================== Plotting ========================

def plot_split(qprod, C, n, rfit_frac=0.9, save_path=None, show=True):
    """
    Plot energy and force decomposition comparing both fitting strategies.
    rfit_frac: fraction of r0 where the rfit parabola is fit (e.g. 0.9).
    """
    r0 = compute_r0(qprod, C, n)
    r_fit = rfit_frac * r0
    E0 = E_full(np.array([r0]), qprod, C, n)[0]
    E_fit = E_full(np.array([r_fit]), qprod, C, n)[0]
    F_fit = F_full(np.array([r_fit]), qprod, C, n)[0]
    k_r0 = n * (n - 1) * C / r0**(n + 2)
    k_rfit = d2E_full(np.array([r_fit]), qprod, C, n)[0]

    # Numerical curvature checks
    dr = 1e-5
    k_r0_num = (E_full(np.array([r0 + dr]), qprod, C, n)[0] - 2 * E0 + E_full(np.array([r0 - dr]), qprod, C, n)[0]) / dr**2
    k_rfit_num = (E_full(np.array([r_fit + dr]), qprod, C, n)[0] - 2 * E_fit + E_full(np.array([r_fit - dr]), qprod, C, n)[0]) / dr**2

    print(f"Parameters: q_i*q_j={qprod:+.3f}  C={C:.3f}  n={n}")
    print(f"  r0            = {r0:.6f}    E(r0) = {E0:.6f}")
    print(f"  r_fit         = {r_fit:.6f}   E(r_fit) = {E_fit:.6f}   F(r_fit) = {F_fit:.6f}")
    print(f"  k at r0       = {k_r0:.6f}  (num: {k_r0_num:.6f})")
    print(f"  k at r_fit    = {k_rfit:.6f}  (num: {k_rfit_num:.6f})")
    print(f"  k ratio       = {k_rfit / k_r0:.2f}x stiffer at r_fit")

    # Range: from ~0.3*r0 to ~3*r0
    r = np.linspace(0.3 * r0, 3.0 * r0, 2000)

    # Compute energies and forces for both strategies
    E       = E_full(r, qprod, C, n)
    E_h_r0  = E_hard_r0(r, qprod, C, n, r0)
    E_s_r0  = E_soft_r0(r, qprod, C, n, r0)
    E_h_rf  = E_hard_rfit(r, qprod, C, n, r_fit)
    E_s_rf  = E_soft_rfit(r, qprod, C, n, r_fit)

    F       = F_full(r, qprod, C, n)
    F_h_r0  = F_hard_r0(r, qprod, C, n, r0)
    F_s_r0  = F_soft_r0(r, qprod, C, n, r0)
    F_h_rf  = F_hard_rfit(r, qprod, C, n, r_fit)
    F_s_rf  = F_soft_rfit(r, qprod, C, n, r_fit)

    fig, (axE, axF) = plt.subplots(2, 1, figsize=(11, 11), sharex=True)

    # --- Energy plot ---
    axE.plot(r, E,            'k-',  linewidth=2.5, label=r'$E(r) = q_i q_j/r + C/r^n$')
    axE.plot(r, E_h_r0 + E0,  'r--', linewidth=1.5, label=r'$E_{\rm hard}^{\rm r0}$ (fit at $r_0$, shifted by $E_{\min}$)')
    axE.plot(r, E_h_rf,       'm-',  linewidth=1.5, label=r'$E_{\rm hard}^{\rm r_{fit}}$ (fit at $0.9 r_0$, with bias)')
    axE.plot(r, E_s_r0,       'b--', linewidth=1.5, label=r'$E_{\rm soft}^{\rm r0} = E - E_{\rm hard}^{\rm r0}$')
    axE.plot(r, E_s_rf,       'g-',  linewidth=1.5, label=r'$E_{\rm soft}^{\rm r_{fit}} = E - E_{\rm hard}^{\rm r_{fit}}$')

    # Mark equilibrium and fitting points
    axE.axvline(r0, color='gray', linestyle=':', alpha=0.7)
    axE.axvline(r_fit, color='purple', linestyle=':', alpha=0.7)
    axE.plot(r0, E0, 'ko', markersize=8, zorder=5)
    axE.plot(r_fit, E_fit, 'mo', markersize=8, zorder=5)
    axE.annotate(f'$r_0$ = {r0:.3f}\n$E_{{\\min}}$ = {E0:.3f}',
                 xy=(r0, E0), xytext=(r0 * 1.25, E0 + 0.15 * abs(E0)),
                 fontsize=10, arrowprops=dict(arrowstyle='->', color='gray'),
                 bbox=dict(boxstyle='round,pad=0.3', fc='lightyellow', ec='gray', alpha=0.8))
    axE.annotate(f'$r_{{fit}}$ = {r_fit:.3f}\n$E$ = {E_fit:.3f}\n$F$ = {F_fit:.3f}',
                 xy=(r_fit, E_fit), xytext=(r_fit * 0.55, E_fit + 0.4 * abs(E0)),
                 fontsize=10, arrowprops=dict(arrowstyle='->', color='purple'),
                 bbox=dict(boxstyle='round,pad=0.3', fc='lavender', ec='purple', alpha=0.8))

    # Shade hard regions
    axE.axvspan(r[0], r_fit, alpha=0.05, color='purple', label=f'hard region ($r < r_{{fit}}={rfit_frac}r_0$)')
    axE.axvspan(r_fit, r0, alpha=0.05, color='orange', label=f'extra soft region ($r_{{fit}} < r < r_0$)')

    axE.set_ylabel('Energy', fontsize=13)
    axE.set_title(f'Hard/Soft Potential Split — r0 fit vs r_fit={rfit_frac}·r0 fit  '
                  f'($q_i q_j$={qprod:+.2f}, $C$={C:.2f}, $n$={n}, '
                  f'$k_{{r0}}$={k_r0:.2f}, $k_{{r_{{fit}}}}$={k_rfit:.2f} = {k_rfit/k_r0:.1f}× stiffer)',
                  fontsize=12)
    axE.legend(loc='upper right', fontsize=9)
    vmin = 1.2 * E0
    vmax = -2.0 * vmin
    axE.set_ylim(vmin, vmax)
    axE.grid(True, alpha=0.3)

    # --- Force plot ---
    axF.plot(r, F,       'k-',  linewidth=2.5, label=r'$F = -dE/dr$')
    axF.plot(r, F_h_r0,  'r--', linewidth=1.5, label=r'$F_{\rm hard}^{\rm r0}$ (fit at $r_0$)')
    axF.plot(r, F_h_rf,  'm-',  linewidth=1.5, label=r'$F_{\rm hard}^{\rm r_{fit}}$ (fit at $0.9 r_0$)')
    axF.plot(r, F_s_r0,  'b--', linewidth=1.5, label=r'$F_{\rm soft}^{\rm r0}$')
    axF.plot(r, F_s_rf,  'g-',  linewidth=1.5, label=r'$F_{\rm soft}^{\rm r_{fit}}$')

    axF.axvline(r0, color='gray', linestyle=':', alpha=0.7)
    axF.axvline(r_fit, color='purple', linestyle=':', alpha=0.7)
    axF.axhline(0, color='gray', linewidth=0.5)
    axF.axvspan(r[0], r_fit, alpha=0.05, color='purple')
    axF.axvspan(r_fit, r0, alpha=0.05, color='orange')

    axF.set_xlabel('$r$ (distance)', fontsize=13)
    axF.set_ylabel('Force', fontsize=13)
    axF.set_title('Force Decomposition', fontsize=14)
    axF.legend(loc='upper right', fontsize=9)
    axF.set_ylim(vmin, vmax)
    axF.grid(True, alpha=0.3)

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Saved figure to {save_path}")

    if show:
        plt.show()

    plt.close(fig)


# ======================== Main ========================

def main():
    parser = argparse.ArgumentParser(description='Plot hard/soft potential split for RigidAtomFF')
    parser.add_argument('--qprod', type=float, default=-1.0,
                        help='Product q_i*q_j (must be < 0 for attractive Coulomb). Default: -1.0')
    parser.add_argument('--C', type=float, default=1.0,
                        help='Repulsion coefficient C. Default: 1.0')
    parser.add_argument('--n', type=int, default=6, choices=[4, 6, 8],
                        help='Repulsion power n. Default: 6')
    parser.add_argument('--rfit', type=float, default=0.9,
                        help='Fraction of r0 where r_fit parabola is fit (e.g. 0.9). Default: 0.9')
    parser.add_argument('--saveFig', type=str, default='',
                        help='Save figure to this path (e.g. hard_soft_split.png)')
    parser.add_argument('--noPlot', action='store_true',
                        help='Do not show the plot (use with --saveFig)')
    args = parser.parse_args()

    if args.qprod >= 0:
        print(f"ERROR: qprod={args.qprod} must be < 0 (attractive Coulomb needed for a minimum)")
        return

    plot_split(
        qprod=args.qprod, C=args.C, n=args.n, rfit_frac=args.rfit,
        save_path=args.saveFig if args.saveFig else None,
        show=not args.noPlot,
    )


if __name__ == '__main__':
    main()
