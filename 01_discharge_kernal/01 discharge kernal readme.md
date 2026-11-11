# 01 — Discharge Kernel (1D Axial Solver)

Solves the 1D axial discharge channel physics using [HallThruster.jl](https://github.com/UM-PEPL/HallThruster.jl).
Produces exit-plane plasma and neutral conditions that feed into the near-field plume solver (Step 2).

**Output at the thruster exit plane (as a function of radius):**
- Ion current density
- Ion energy distribution
- Neutral density

---

## Folder Structure

_NOTE: Serial numbers next to the codes and documentation indicate the sequence of simulation._
```
01_discharge_kernal/
├── Code - Julia Suite/
│   ├── 1_Code_Julia_Suite_README.md         ← usage and run order
│   ├── 2_My_HET.jl                          ← core simulation script
│   └── 3_2D axisymmetric cross section.jl   ← data extraction for visualiser
└── Code documentation/
    ├── 1_Discharge Kernel README.md         ← this file
    ├── 2_myHET_7mN.md                       ← My_HET.jl walkthrough
    └── 3_2D_axisymmetric_cross_section.md   ← extraction script walkthrough
```

---

## Code Documentation

| File | Covers | Description |
|---|---|---|
| `myHET_7mN.md` | `My_HET.jl` | Full walkthrough — geometry, B-field, config, thrust sweep, plots, Excel export, design scaling reference |
| `2D_axisymmetric_cross_section.md` | `2D axisymmetric cross section.jl` | Data extraction walkthrough — plasma profiles, viz_data bundle, operating conditions |

---

## Dependencies

- [HallThruster.jl](https://github.com/UM-PEPL/HallThruster.jl) v0.23.2
- Julia 1.12.6
- Plots.jl, XLSX.jl, Printf, Statistics
