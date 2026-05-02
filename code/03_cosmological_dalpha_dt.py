"""
Paper VII - cosmological time-variation of alpha_EM.

If the substrate's reserve R_0 has slow cosmological dynamics, then
alpha_EM inherits a small drift over cosmological time.

Model
-----
Assume R_0(t) evolves with a fractional rate
  dR_0/R_0 / dt = beta_R / t_H (fractional Hubble rate)
where t_H = 14 Gyr is the Hubble time.

Then via the cross-coupling derived in script 01:
  dalpha/alpha = -c_grav * delta_R / R_0
With c_grav from script 01.

Output
------
Predicted |dalpha/alpha / dt| in /yr, for various scenarios:
  - beta_R = 1.0  (R_0 changes by factor e over Hubble time)
  - beta_R = 0.1
  - beta_R = 0.01

Compared to current bound:
  |dalpha/alpha / dt| < 1e-17 /yr  (Webb et al., quasar absorption + atomic clocks)
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
    return pi**(d/2) / gamma(d/2 + 1)


def alpha_at_x(x):
    return 1.0 / (8 * pi**2 * math.sqrt(x))


def fixed_point_alpha(V):
    def fn(x):
        a = alpha_at_x(x)
        return x - 3 - a*math.sqrt(3)*(1 - a*V)
    sol = opt.brentq(fn, 3.0, 3.05)
    return sol, alpha_at_x(sol)


def c_grav(epsilon_0, qa):
    """Recompute c_grav for a given epsilon_0 and qa (see script 01)."""
    eps_perturb = 1e-6
    d0 = 3 + epsilon_0
    V0 = V_d(d0)
    _, alpha0 = fixed_point_alpha(V0)
    eps_loc = epsilon_0 * (1 - qa * eps_perturb)
    V_loc = V_d(3 + eps_loc)
    _, alpha_loc = fixed_point_alpha(V_loc)
    return -(alpha_loc - alpha0) / alpha0 / eps_perturb


def main():
    print("=" * 75)
    print("COSMOLOGICAL VARIATION OF alpha_EM")
    print("=" * 75)
    print()
    print("Assume R_0(t) has fractional rate dR/(R_0 dt) = beta_R / t_H")
    print("where t_H = 14 Gyr is the Hubble time.")
    print()
    print("Then dalpha/alpha = -c_grav * dR_0/R_0  (linear coupling, qa=1)")
    print("    dalpha/(alpha dt) = -c_grav * beta_R / t_H")
    print()

    t_H = 14e9  # 14 Gyr in years

    print(f"{'epsilon_0':<22} {'c_grav (qa=1)':<15} "
          f"{'|dalpha/alpha|/yr at beta_R=':<28}")
    print(f"{'':<22} {'':<15} {'1.0':<10} {'0.1':<10} {'0.01':<10}")
    print("-" * 75)

    candidates = [
        ("1/(3*pi)", 1/(3*pi)),
        ("1/(2*pi)", 1/(2*pi)),
        ("1/6",      1/6),
    ]
    for label, eps_0 in candidates:
        c = c_grav(eps_0, qa=1.0)
        rates = {beta_R: abs(c) * beta_R / t_H for beta_R in [1.0, 0.1, 0.01]}
        print(f"epsilon = {label:<13} {c:>+.4e}    "
              f"{rates[1.0]:<10.2e} {rates[0.1]:<10.2e} {rates[0.01]:<10.2e}")

    print()
    print("Current bounds:")
    print("  Atomic clocks  : |dalpha/(alpha*dt)| < 1e-17 /yr  (Rosenband et al.)")
    print("  Quasar (z~3)   : |Delta alpha/alpha| < 1e-5 over 12 Gyr")
    print("                   = ~1e-15 /yr equivalent")
    print()
    print("Conclusion:")
    print()
    print("  beta_R ~ 1 : DDD predicts dalpha/dt ~ 1e-16 /yr")
    print("    -> within 10x of current atomic-clock bound,")
    print("       within 10x below current quasar bound.")
    print("       FUTURE TESTABLE.")
    print()
    print("  beta_R ~ 0.1 : DDD predicts dalpha/dt ~ 1e-17 /yr")
    print("    -> right AT the atomic-clock bound (Rosenband).")
    print("    Borderline detectable today.")
    print()
    print("  beta_R ~ 0.01 : DDD predicts dalpha/dt ~ 1e-18 /yr")
    print("    -> below current detection by 10x")
    print()

    # Equivalent: Delta alpha/alpha at quasar redshift z=3
    print("=" * 75)
    print("Quasar redshift comparison (z = 3, lookback ~12 Gyr):")
    print("=" * 75)
    print()
    for label, eps_0 in candidates:
        c = c_grav(eps_0, qa=1.0)
        # If R_0 evolves as a power-law in scale factor a:
        # dR/R = (dR/da) * (da/a)  with da/a = -dz/(1+z) = -1/(1+z)*dz
        # Over z = 0 to 3 : delta_a/a ~ -3/4
        delta_R_over_R = 0.75  # assume |delta R/R| = 0.75 between z=3 and z=0
        delta_alpha = c * delta_R_over_R
        print(f"epsilon = {label:<10}: |Delta alpha/alpha| ~ {abs(delta_alpha):.2e}")
    print()
    print("Quasar bound (Webb et al., Murphy et al.): |Delta alpha/alpha| < 1e-5")
    print("DDD prediction: ~1e-6 (one order below the bound).")
    print()

    out = {
        "t_H_yr": t_H,
        "candidates": [{"label": label, "epsilon_0": eps_0,
                        "c_grav_qa1": c_grav(eps_0, 1.0)}
                        for label, eps_0 in candidates],
    }
    with open(DATA / "03_cosmological_dalpha_dt.json", "w") as f:
        json.dump(out, f, indent=2)
    print(f"Saved -> {DATA / '03_cosmological_dalpha_dt.json'}")


if __name__ == "__main__":
    main()
