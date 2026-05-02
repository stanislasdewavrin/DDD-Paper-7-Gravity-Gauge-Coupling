"""
Paper VII - analytical derivation of c_grav prefactor.

Setup
-----
Paper II : near a mass M, R(r)/R_0 = sqrt(1 - beta/r)
           with beta = 2GM/c^2

Paper VI : alpha_EM = 1/(8*pi^2*sqrt(x*)) ;  x* = 3 + alpha*sqrt(3)*(1 - alpha*V)
           V = pi^(d/2)/Gamma(d/2+1), d = 3 + epsilon
           epsilon comes from small-world correction whose density scales with R

Cross-coupling
--------------
Hypothesis: shortcut density rho_short(R) ~ R^q for small q.
=> epsilon_local(R) = epsilon_0 * (R/R_0)^(q*a)  where a is small-world exponent.

For R^2/R_0^2 = 1 - beta/r, we have R/R_0 ~ 1 - beta/(2r) ~ 1 - |Phi|/c^2.

So epsilon_local ~ epsilon_0 * (1 - q*a*|Phi|/c^2)
   d_eff_local  = 3 + epsilon_local
   V_local      = V_0 + (dV/dd)|_d_0 * delta_d

Solving the fixed-point self-consistently with V_local gives
delta_alpha/alpha = -c_grav * |Phi|/c^2.

This script computes c_grav for various q*a values.
"""

import math
from math import pi, gamma
import scipy.optimize as opt
from scipy.special import digamma
import numpy as np
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
DATA = HERE / "data"; DATA.mkdir(exist_ok=True)


def V_d(d):
    """Volume of unit ball in dimension d."""
    return pi**(d/2) / gamma(d/2 + 1)


def dV_dd(d):
    """Derivative of V_d w.r.t. d."""
    # V(d) = pi^(d/2) / Gamma(d/2 + 1)
    # ln V = (d/2) ln pi - ln Gamma(d/2+1)
    # d(ln V)/dd = (1/2) ln pi - (1/2) psi(d/2+1)
    # dV/dd = V * (1/2) * (ln pi - psi(d/2+1))
    return V_d(d) * 0.5 * (math.log(pi) - digamma(d/2 + 1))


def alpha_at_x(x):
    return 1.0 / (8 * pi**2 * math.sqrt(x))


def fixed_point_alpha(V):
    """Solve self-consistent equation for alpha_EM given V."""
    def fn(x):
        a = alpha_at_x(x)
        return x - 3 - a*math.sqrt(3)*(1 - a*V)
    sol = opt.brentq(fn, 3.0, 3.05)
    return sol, alpha_at_x(sol)


def c_grav(epsilon_0, qa):
    """
    Compute c_grav such that delta_alpha/alpha = -c_grav * |Phi|/c^2.

    Procedure:
      1. Compute baseline alpha_0 with d_eff_0 = 3 + epsilon_0
      2. Compute alpha at d_eff = 3 + epsilon_0 * (1 - qa * eps_perturb)
         for small eps_perturb representing |Phi|/c^2
      3. delta_alpha/alpha = (alpha_perturbed - alpha_0) / alpha_0
      4. c_grav = -(delta_alpha/alpha) / eps_perturb
    """
    d0 = 3 + epsilon_0
    V0 = V_d(d0)
    x0, alpha0 = fixed_point_alpha(V0)

    eps_perturb = 1e-6
    epsilon_local = epsilon_0 * (1 - qa * eps_perturb)
    d_local = 3 + epsilon_local
    V_local = V_d(d_local)
    x_local, alpha_local = fixed_point_alpha(V_local)

    delta_alpha = alpha_local - alpha0
    delta_alpha_over_alpha = delta_alpha / alpha0
    c = -delta_alpha_over_alpha / eps_perturb
    return c, alpha0, V0


def main():
    print("=" * 70)
    print("DERIVATION OF c_grav")
    print("=" * 70)
    print()
    print("Hypothesis: epsilon_local(R) = epsilon_0 * (R/R_0)^(q*a)")
    print("            with R/R_0 = 1 - |Phi|/c^2 (linear in Phi)")
    print()
    print("Computing c_grav for three epsilon_0 candidates and various qa:")
    print()

    candidates = [
        ("eps = 1/(3*pi)", 1/(3*pi)),
        ("eps = 1/(2*pi)", 1/(2*pi)),
        ("eps = 1/6",      1/6),
    ]

    for label, eps_0 in candidates:
        d0 = 3 + eps_0
        V0 = V_d(d0)
        dVdd0 = dV_dd(d0)
        print(f"--- {label}, d_eff_0 = {d0:.5f}, V = {V0:.5f}, dV/dd = {dVdd0:.5f} ---")
        for qa in [0.5, 1.0, 2.0, 4.0]:
            c, alpha0, V0 = c_grav(eps_0, qa)
            # Predict delta_alpha/alpha for various physical situations
            phi_earth = 6.95e-10  # |Phi|/c^2 at Earth surface
            phi_sun = 2.12e-6     # |Phi|/c^2 at Sun surface
            phi_ns = 0.1          # |Phi|/c^2 at neutron star surface
            print(f"  qa={qa:>4}: c_grav = {c:+.4e}  "
                  f"=> delta_alpha/alpha (Earth) = {-c*phi_earth:+.2e}, "
                  f"(Sun) = {-c*phi_sun:+.2e}, (NS) = {-c*phi_ns:+.2e}")
        print()

    # Most likely scenario: linear coupling (qa = 1)
    print("=" * 70)
    print("Linear coupling (qa = 1): the natural baseline")
    print("=" * 70)
    print()
    qa = 1.0
    print(f"Prediction: delta_alpha/alpha = -c_grav * |Phi|/c^2")
    print()
    for label, eps_0 in candidates:
        c, _, _ = c_grav(eps_0, qa)
        # Atomic clock altitude test (1 km elevation difference)
        delta_phi_altitude = 9.81 * 1000 / (3e8)**2
        delta_alpha_alt = abs(-c * delta_phi_altitude)
        # Atomic clock at 1 AU vs Earth surface (annual variation due to eccentricity)
        delta_phi_orbit = 5e-10  # ~ 0.017 * (GM_sun/r_AU/c^2)
        delta_alpha_orbit = abs(-c * delta_phi_orbit)
        # Quasar bound on cosmological variation: |delta_alpha/alpha| < 1e-5
        # over redshift z ~ 1, where Phi structure ~ GM_universe/r_obs ~ 1e-5
        print(f"{label}: c_grav = {c:.4e}")
        print(f"  Atomic clock altitude (h=1km) test : |delta_alpha/alpha| = {delta_alpha_alt:.2e}")
        print(f"  Atomic clock orbital variation    : |delta_alpha/alpha| = {delta_alpha_orbit:.2e}")
        print()

    print("Current bounds:")
    print("  Atomic clocks (Cs vs Sr): |dalpha/alpha| < 1e-17 over 1 year")
    print("  Quasar absorption (z~3) : |Delta alpha/alpha| < 1e-5")
    print("  GPS satellite vs ground : |delta_alpha/alpha| < 1e-15")
    print()

    # Save results
    out = {}
    for label, eps_0 in candidates:
        for qa in [0.5, 1.0, 2.0]:
            c, _, _ = c_grav(eps_0, qa)
            out[f"{label}, qa={qa}"] = c
    with open(DATA / "01_cgrav.json", "w") as f:
        json.dump(out, f, indent=2)
    print(f"Saved -> {DATA / '01_cgrav.json'}")


if __name__ == "__main__":
    main()
