# 01 — Discharge Kernel (1D Axial Solver)

Solves the 1D axial discharge channel physics using [HallThruster.jl](https://github.com/UM-PEPL/HallThruster.jl).
Produces exit-plane plasma and neutral conditions that feed into the near-field plume solver (Step 2).

**Outputs at the thruster exit plane (as a function of radius):**
- Ion current density
- Ion energy distribution
- Neutral density

---

## Folder Structure

```
01_discharge_kernal/
├── Code - Julia Suite/
│   └── My_HET.jl           ← simulation script
└── Code documentation/
    ├── Discharge Kernal readme.md   ← module walkthrough
    └── myHET_7mN.md                ← full code documentation
```

---

## What This Module Does

This is **Step 1** of the five-step HET plume–spacecraft interaction framework.
The 1D discharge solver computes the internal channel plasma physics — ionisation, electron heating, and ion acceleration — and delivers the exit-plane conditions needed for the downstream plume expansion model.

---

## Dependencies

- [HallThruster.jl](https://github.com/UM-PEPL/HallThruster.jl) v0.23.2
- Julia 1.12.6
- Plots.jl, XLSX.jl, Printf, Statistics
