# Code — Julia Suite

Julia simulation scripts for the 1D Hall Effect Thruster discharge kernel.
Run these in order in a single REPL session — each script reads variables left in scope by the previous one.

---

## Files

| File | Description |
|---|---|
| `My_HET.jl` | Core simulation — geometry, B-field, Config, SimParams, baseline run, thrust sweep, 4-panel plots, Excel export |

---

## Usage

```julia
# Run from the Julia REPL
julia> include("My_HET.jl")
```

Expected console output:

```
Running simulation...
── Results ───────────────────
Thrust            : 7.xx mN
Discharge current : x.xx A
Power             : xxx.x W
Isp               : xxxx s
Anode efficiency  : xx.x %
```

## Outputs

| File | Description |
|---|---|
| `MyHET_7mN_results.png` | 4-panel plot: B-field, electron temperature, plasma potential, ion current density |
| `MyHET_7mN_results.xlsx` | Plasma properties · Performance summary · Thrust sweep · Magnetic field profile |

---

## Dependencies

```julia
using Pkg
Pkg.add(["HallThruster", "Plots", "XLSX"])
```

| Package | Version | Purpose |
|---|---|---|
| HallThruster.jl | 0.23.2 | 1D Hall thruster discharge solver |
| Plots.jl | — | Visualisation (GR backend) |
| XLSX.jl | — | Excel export |
| Printf | stdlib | Formatted console output |
| Statistics | stdlib | Time-averaging |
