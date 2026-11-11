# 01 — Discharge Kernel (1D Axial Solver)

Solves the 1D axial discharge channel physics using [HallThruster.jl](https://github.com/UM-PEPL/HallThruster.jl).
Produces exit-plane plasma and neutral conditions that feed into the near-field plume solver (Step 2).

**Output at the thruster exit plane (as a function of radius):**
- Ion current density
- Ion energy distribution
- Neutral density

---

## Folder Structure

```
01_discharge_kernal/
├── Code - Julia Suite/
│   ├── 1_Code_Julia_Suite_README.md       ← usage and run order
│   ├── My_HET.jl                          ← core simulation script
│   └── 2D axisymmetric cross section.jl   ← data extraction for visualiser
└── Code documentation/
    ├── 1 Discharge Kernel README.md       ← this file
    ├── myHET_7mN.md                       ← My_HET.jl walkthrough
    └── 2D_axisymmetric_cross_section.md   ← extraction script walkthrough
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
