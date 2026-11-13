# 1 — Discharge Kernel (1D Axial Solver)

Solves the 1D axial discharge channel physics using [HallThruster.jl](https://github.com/UM-PEPL/HallThruster.jl).
Produces exit-plane plasma and neutral conditions that feed into the near-field plume solver (Step 2).

**Output at the thruster exit plane (as a function of radius):**
- Ion current density
- Ion energy distribution
- Neutral density

---

## Folder Structure

> **Note:** Serial numbers on files indicate the sequence in which scripts should be run and documentation should be read.

```
01_discharge_kernal/
├── Code - Julia Suite/
|   ├── 0_ConfigV3.jl                             ← Automated config generator
│   ├── 1_My_HET.jl                               ← core simulation script
│   ├── 2_2D axisymmetric cross section.jl        ← data extraction for visualiser
│   ├── 3_Visualise thruster 2d.jl                ← 2D axisymmetric visualiser
│   ├── 4_Derived_Quantities.jl                   ← derived plasma physics + Excel export
|   └── README.md                                 ← usage and run order
└── Code documentation/
    ├── 0_ConfigV3.md                             ← Automated config generator
    ├── 1_myHET_7mN.md                            ← My_HET.jl walkthrough
    ├── 2_2D_axisymmetric_cross_section.md        ← extraction script walkthrough
    ├── 3_Visualise_thruster_2d.md                ← 2D visualiser walkthrough
    ├── 4_Derived_Quantities.md                   ← derived quantities walkthrough
    └── README.md                                 ← this file
```

---

## Code Documentation

| # | File | Covers | Description |
|---|---|---|---|
| 0 | `0_ConfigV3.md` | `0_ConfigV3.jl` | Automated CONFIG generator — Nelder-Mead + bisection search for optimal (Vd, mdot, geometry). Run standalone before the main pipeline. |
| 1 | `1_myHET_7mN.md` | `2_My_HET.jl` | Full walkthrough — geometry, B-field, config, thrust sweep, plots, Excel export, design scaling reference |
| 2 | `2_2D_axisymmetric_cross_section.md` | `3_2D axisymmetric cross section.jl` | Data extraction — plasma profiles, viz_data bundle, operating conditions |
| 3 | `3_Visualise_thruster_2d.md` | `4_Visualise thruster 2d.jl` | 8-panel 2D visualiser — zone shading, plume divergence cones, geometry schematic |
| 4 | `4_Derived_Quantities.md` | `5_Derived_Quantities.jl` | Derived plasma physics — E-field, ion velocity, Hall parameter, E×B drift, Excel export |

---

## Dependencies

- [HallThruster.jl](https://github.com/UM-PEPL/HallThruster.jl) v0.23.2
- Julia 1.12.6
- Plots.jl, XLSX.jl, Printf, Statistics
