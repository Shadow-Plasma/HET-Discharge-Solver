# 01 — Discharge Kernel (1D Axial Solver)

Solves the 1D axial discharge channel physics using [HallThruster.jl](https://github.com/UM-PEPL/HallThruster.jl).
Produces exit-plane plasma and neutral conditions that feed into the near-field plume solver (Step 2).

**Output at the thruster exit plane (as a function of radius):**
- Ion current density
- Ion energy distribution
- Neutral density

---

## Files

| File | Description |
|---|---|
| `MyHET_7mN.jl` | Simulation script — geometry, B-field, config, thrust sweep, plots, Excel export |
| `MyHET_7mN.md` | Full code walkthrough and design scaling reference |

---

## Dependencies

- [HallThruster.jl](https://github.com/UM-PEPL/HallThruster.jl) v0.23.2
- Julia 1.12.6
