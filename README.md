# Paper VII — Gravity-Gauge Coupling and the Joint Wilson-Loop Back-Reaction

This paper proposes that gravity (Paper II) and electromagnetism (Paper VI)
are two back-reaction channels of the same discrete substrate. On the unit
cube, two distinct closed Wilson objects coexist:

- **Body-diagonal Wilson loop** (squared length 3) → α_EM
- **Cube-surface Wilson surface** (area 6) → α_G at the substrate scale

Under the assumption that both sectors share the same combinatorial prefactor
8π² and one-loop structure, the leading structural ratio is

```
α_G^lat / α_EM^lat = √(3/6) = 1/√2 ≈ 0.707
```

## The mystery this paper addresses

Why is gravity ~40 orders of magnitude weaker than electromagnetism for
ordinary particles? Paper VII reframes this as a **scale hierarchy**:
at the substrate scale, the two couplings are geometrically comparable
(within a factor √2). The observed weakness comes from the fact that
ordinary particle masses lie far below the substrate mass scale.

## Headline results

- **Tree-level structural ratio** α_G/α_EM = 1/√2 (conditional on shared
  Wilson-loop prefactor)
- **Phenomenological substrate mass scale** m_lat ≈ 0.072 M_P,
  ℓ_lat ≈ 14 ℓ_P (using α_EM(0) as low-energy anchor;
  full RG treatment deferred to Paper XV)
- **Cross-coupling** c_grav ≈ 1.4–2.2 × 10⁻⁶ between the local
  reserve R(r) and α_EM, under the minimal linear ansatz
  ε_local(R) = ε_0 R/R_0
- **Cosmological signature** |dα/α / dt| ≈ 10⁻¹⁶/yr at β_R ∼ 1, in
  tension with current Rosenband bound 10⁻¹⁷/yr unless β_R ≲ 0.1
  or weaker reserve dependence

## Falsifiers

- Atomic clock comparisons across gravitational potentials (within
  factor ~10 of present sensitivity, within minimal linear ansatz)
- Cosmological drift of α_EM near current atomic-clock bounds
- Quasar absorption Δα/α ~ 10⁻⁶ at z=3, below current bound 10⁻⁵

## Repository structure

```
paper.tex          — main LaTeX source
paper.pdf          — compiled paper
references.bib     — bibliography (Rosenband, Webb, Murphy, Wiens-Tuller)
code/              — Python scripts:
   01_cgrav_analytical.py         — derives c_grav by self-consistent
                                     fixed-point with dV/dd via digamma
   02_cross_coupling_simulation.py — lattice L=40 with central mass,
                                     measures BFS dimension shift
   03_cosmological_dalpha_dt.py   — projects cosmological drift vs
                                     Rosenband and Webb bounds
data/              — JSON output from the scripts
figures/           — figures used in the paper
```

## Reproducing the analytics

```bash
cd code
python 01_cgrav_analytical.py        # → data/01_cgrav.json
python 02_cross_coupling_simulation.py  # → data/02_cross_coupling.json
python 03_cosmological_dalpha_dt.py  # → data/03_cosmological_dalpha_dt.json
```

## Compiling the paper

```bash
pdflatex paper.tex
bibtex paper
pdflatex paper.tex
pdflatex paper.tex
```

## Companion papers

- **Paper II** (Schwarzschild from drainage back-reaction):
  https://github.com/stanislasdewavrin/DDD-Paper-2-gravity
- **Paper VI** (α_EM fixed-point formula):
  https://github.com/stanislasdewavrin/DDD-Paper-6-electromagnetism

## Citation

Stanislas Dewavrin, *Discrete Drainage Dynamics VII: Gravity-Gauge
Coupling and the Joint Wilson-Loop Back-Reaction*, May 2026.
