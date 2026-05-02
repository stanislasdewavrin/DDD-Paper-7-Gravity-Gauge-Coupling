"""
Paper VII - numerical simulation of gravity-gauge cross-coupling.

Setup
-----
Lattice with central mass:
  - 3D cubic L^3 lattice
  - Local reserve R(r) = R_0 * sqrt(1 - beta/r)  (Paper II Schwarzschild)
  - Body-diagonal shortcut density: rho_local(r) ~ R(r)^q * mu_0
    (depleted near central mass)

Procedure
---------
For each beta value (mass strength):
  1. Compute R(r) field on lattice
  2. Generate Poisson(mu_local) shortcuts with mu_local depending on R(r)
  3. Measure BFS expansion dimension d_BFS in two regions:
     (a) far from center (control)
     (b) close to center (depleted region)
  4. Extract delta_d_BFS / d_BFS

If the cross-coupling hypothesis is correct, we should see:
  delta_d_BFS proportional to <|Phi|/c^2> in the depleted region.
"""
import numpy as np
import json
from collections import deque
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
DATA = HERE / "data"; DATA.mkdir(exist_ok=True)

DIRS4 = np.array([
    [+1, +1, +1],
    [+1, +1, -1],
    [+1, -1, +1],
    [-1, +1, +1],
], dtype=np.int32)


def idx_of(i, j, k, L):
    return ((i % L) * L + (j % L)) * L + (k % L)


def build_substrate_with_mass(L, beta=0.0, mu_0=1.0, q=1.0, seed=42):
    """
    Build cubic lattice + body-diagonal Poisson shortcuts with local density
    modulated by R(r) = R_0 * sqrt(1 - beta/r), where r = distance from center.

    Returns
    -------
    adj : list of lists (adjacency)
    R_field : array (L,L,L) of local reserves
    """
    rng = np.random.default_rng(seed)
    N = L**3
    cx = cy = cz = L // 2

    # Distance from center for each node
    R_field = np.ones((L, L, L))
    for i in range(L):
        for j in range(L):
            for k in range(L):
                r = np.sqrt((i-cx)**2 + (j-cy)**2 + (k-cz)**2) + 1e-6
                ratio_sq = max(0.01, 1 - beta/r)  # avoid going negative
                R_field[i, j, k] = np.sqrt(ratio_sq)

    # Build cubic adjacency
    adj = [[] for _ in range(N)]
    for i in range(L):
        for j in range(L):
            for k in range(L):
                u = idx_of(i, j, k, L)
                for di, dj, dk in [(1, 0, 0), (0, 1, 0), (0, 0, 1)]:
                    v = idx_of(i+di, j+dj, k+dk, L)
                    adj[u].append(v); adj[v].append(u)

    # Poisson body-diagonal shortcuts with density modulated by R(r)^q
    for i in range(L):
        for j in range(L):
            for k in range(L):
                u = idx_of(i, j, k, L)
                R_local = R_field[i, j, k]
                mu_local = mu_0 * R_local**q
                k_short = rng.poisson(mu_local)
                for _ in range(k_short):
                    didx = int(rng.integers(0, 4))
                    dvec = DIRS4[didx]
                    sign = 1 if rng.random() < 0.5 else -1
                    di, dj, dk = (sign*dvec[0], sign*dvec[1], sign*dvec[2])
                    v = idx_of(i+di, j+dj, k+dk, L)
                    if v != u:
                        adj[u].append(v); adj[v].append(u)
    return adj, R_field


def bfs_dim(adj, source, r_min=2, r_frac_max=0.4):
    N = len(adj)
    dist = [-1] * N
    dist[source] = 0
    q = deque([source])
    while q:
        u = q.popleft()
        for v in adj[u]:
            if dist[v] < 0:
                dist[v] = dist[u] + 1
                q.append(v)
    max_d = max(dist)
    N_at = np.zeros(max_d + 1, dtype=int)
    for d in dist:
        if d >= 0:
            N_at[d] += 1
    cumul = np.cumsum(N_at)
    rs = np.arange(len(cumul))
    r_max = int(r_frac_max * len(cumul))
    mask = (rs >= r_min) & (rs <= r_max) & (cumul > 0)
    log_r = np.log(rs[mask])
    log_N = np.log(cumul[mask])
    d_H, _ = np.polyfit(log_r, log_N, 1)
    return float(d_H)


def main():
    L = 40
    print(f"L = {L}")
    print()

    rows = []
    # Sweep beta values: from 0 (no mass) to L/8 (significant gravitational depletion)
    betas = [0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0]
    print(f"{'beta':<8} {'<R/R_0> at center':<18} {'<|Phi|/c^2> proxy':<18} {'d_BFS center':<14} {'d_BFS edge':<14} {'shift':<12}")
    print("-" * 95)

    for beta in betas:
        # Use 3 seeds for stats
        d_BFS_centers = []
        d_BFS_edges = []
        for seed in range(3):
            adj, R_field = build_substrate_with_mass(L, beta=beta, seed=seed)
            cx = L // 2
            # BFS from center
            src_center = idx_of(cx, cx, cx, L)
            d_BFS_centers.append(bfs_dim(adj, src_center))
            # BFS from edge (corner)
            src_edge = idx_of(2, 2, 2, L)
            d_BFS_edges.append(bfs_dim(adj, src_edge))

        # Average <R/R_0> at center
        R_at_center = R_field[cx-3:cx+4, cx-3:cx+4, cx-3:cx+4].mean()
        # |Phi|/c^2 proxy
        phi_proxy = beta / max(2, 5.0)  # approximate at r ~ 5

        d_c = np.mean(d_BFS_centers)
        d_e = np.mean(d_BFS_edges)
        shift = d_c - d_e

        print(f"{beta:<8.2f} {R_at_center:<18.5f} {phi_proxy:<18.5f} "
              f"{d_c:<14.5f} {d_e:<14.5f} {shift:<+12.5f}")
        rows.append({
            "beta": beta,
            "R_at_center": float(R_at_center),
            "phi_proxy": float(phi_proxy),
            "d_BFS_center": d_c,
            "d_BFS_edge": d_e,
            "shift": shift,
        })

    print()
    print("Interpretation:")
    print("  - At beta = 0: no mass, no depletion, both d_BFS values should be similar")
    print("  - At beta > 0: depletion at center reduces shortcut density,")
    print("    so d_BFS_center < d_BFS_edge (shift negative)")
    print()
    print(f"Linear fit shift vs beta: ", end="")
    bs = np.array([r["beta"] for r in rows])
    ss = np.array([r["shift"] for r in rows])
    slope, intercept = np.polyfit(bs, ss, 1)
    print(f"slope = {slope:+.5f}, intercept = {intercept:+.5f}")

    out = {"L": L, "rows": rows, "linear_fit": {"slope": slope, "intercept": intercept}}
    with open(DATA / "02_cross_coupling.json", "w") as f:
        json.dump(out, f, indent=2)


if __name__ == "__main__":
    main()
