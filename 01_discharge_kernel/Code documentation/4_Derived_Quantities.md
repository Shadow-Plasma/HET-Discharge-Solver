# Derived Quantities — Plasma Physics Analysis

Computes all derived plasma quantities from `viz_data` — no new simulation is run.
Produces a 6-panel derived-quantities figure and an 8-sheet Excel workbook with live formula columns.

> **Prerequisites:** Run `My_HET.jl` and `2D axisymmetric cross section.jl` first in the same REPL session. This script reads `viz_data` from memory — it will fail in a fresh REPL.

---

## Usage

```julia
# Run all scripts in order in the same REPL session
julia> include("My_HET.jl")
julia> include("2D axisymmetric cross section.jl")
julia> include("Visualise thruster 2d.jl")
julia> include("Derived_Quantities.jl")
```

**Output files:**
- `MyHET_derived.png` — 6-panel derived plasma physics plot
- `MyHET_derived_data.xlsx` — 8-sheet Excel workbook with live formulas

---

## Physical Constants

| Constant | Symbol | Value | Unit |
|---|---|---|---|
| Elementary charge | `e_q` | 1.60218 × 10⁻¹⁹ | C |
| Electron mass | `me` | 9.10938 × 10⁻³¹ | kg |
| Xe atomic mass | `Mi` | 2.18014 × 10⁻²⁵ | kg |
| Standard gravity | `g0` | 9.80665 | m/s² |
| e-n transfer coefficient | `k_en` | 2.5 × 10⁻¹³ | m³/s |

> `Mi` uses the natural Xe standard atomic weight (131.293 u). Xe-131 exact isotope = 130.905 u → 2.1737 × 10⁻²⁵ kg.

---

## Derived Quantities

### 1. Electric Field E(z)

Computed by central-difference differentiation of the plasma potential:

```
E(z) = -dφ/dz   [V/m]
```

Edge rows use forward/backward differences. Positive E points toward the anode.

---

### 2. Ion Velocity v_ion(z)

Derived from energy conservation — ions accelerated from the potential maximum to position z:

```
v_ion(z) = sqrt( 2 · e_q · (φ_max − φ(z)) / Mi )   [m/s]
```

Negative values of `(φ_max − φ(z))` are clipped to zero to suppress numerical noise near the anode.

---

### 3. Ionisation Rate Coefficient k_iz(Te)

Piecewise fit from Goebel & Katz (2008), Appendix E — Table E-1. Valid for Xe, Te ≈ 0.5–10 eV. Accuracy ~2–5% vs tabulated values.

```
v̄ₑ = sqrt( 8 · e_q · Te / (π · me) )      [Maxwellian mean speed]

Te < 5 eV:   k_iz = 1e-20 · (3.97 + 0.643·Te − 0.0368·Te²) · exp(−12.127/Te) · v̄ₑ
Te ≥ 5 eV:   k_iz = 1e-20 · (−1.031e-4·Te² + 6.386·exp(−12.127/Te)) · v̄ₑ
```

Te is floored at 0.5 eV (formula lower bound).

---

### 4. Ionisation Source Term S_iz(z)

```
S_iz(z) = k_iz(Te) · nn(z) · ni(z)   [m⁻³ s⁻¹]
```

Quasineutrality assumed: `ne ≈ ni`.

---

### 5. Hall Parameter Ω_e(z)

```
Ω_e(z) = ω_ce / ν_en

ω_ce = e_q · B / me          [electron cyclotron frequency]
ν_en  = k_en · nn            [electron-neutral collision frequency]
```

`ν_en` is interpolated from the plasma cell grid onto the B-field grid. A floor of +1 is added to the denominator to prevent division by zero in the plume.

> **Note:** This is an order-of-magnitude estimate. Exact values require measured e-n cross-section data.

---

### 6. Electron Larmor Radius r_L(z)

```
r_L(z) = sqrt( 2 · me · Te_J ) / ( e_q · B + 1e-10 )   [m]
```

Te is interpolated onto the B-field grid. The 1e-10 floor prevents division by zero in the plume where B → 0.

`r_L << channel_width` confirms electrons are magnetised throughout the channel.

---

### 7. E×B Drift Velocity v_ExB(z)

```
v_ExB(z) = |E(z)| / B(z)   [m/s]
```

B is interpolated from the B-field grid onto the plasma cell grid.

---

## Scalar Performance Metrics

Printed to the console after the derivations:

| Metric | Expression |
|---|---|
| Acceleration potential Δφ | `φ_max − φ_exit` |
| S_iz-weighted exit velocity | `Σ(weight_i · v_ion_i)` where `weight_i = S_iz_i / ΣS_iz` |
| Isp (S_iz-weighted) | `v_exit / g0` |
| Thrust estimate | `mdot_ion · v_exit` |
| Propellant utilisation η_m | `mdot_ion / mdot` |
| Beam voltage fraction η_v | `Δφ / Vd` |
| Beam efficiency η_beam | `η_m × η_v` |
| Peak E-field | `max(|E(z)|)` and location |
| S_iz peak and FWHM | Location and width of ionisation zone |
| Peak Hall parameter | `max(Ω_e)` |

---

## Panels

| Panel | Quantity | Colour |
|---|---|---|
| 1 | Electric field `E(z) = −dφ/dz` | Yellow |
| 2 | Ion velocity `v_ion(z)` from energy conservation | Purple |
| 3 | Ionisation source `S_iz(z) = k_iz · nn · ni` | Green |
| 4 | Hall parameter `Ω_e(z)` | Pink |
| 5 | Electron Larmor radius `r_L(z)` | Light blue |
| 6 | E×B drift velocity `v_ExB(z)` | Cyan |

All panels carry zone shading (ionisation: yellow, acceleration: blue) and an exit-plane dashed line.

---

## Excel Export — 8-Sheet Workbook

Design principle: each data sheet holds only **primitive simulation outputs**. Every column that can be derived from primitives is an **Excel formula** — changing a primitive cell propagates through the whole workbook automatically.

| Sheet | Contents |
|---|---|
| `Constants` | Physical constants and operating conditions — referenced by all formula sheets |
| `Electric_Field` | z, φ(z) primitives; E(z) = −dφ/dz formula column; Julia cross-check column |
| `Ion_Velocity` | z, φ(z) primitives; φ_max, dV, v_ion formula columns |
| `Ionisation_Source` | z, Te, nn, ni primitives; k_iz (Goebel & Katz), S_iz, weight formula columns |
| `Hall_Parameter` | z_B, B, ν_en primitives; ω_ce, Ω_e formula columns |
| `Larmor_Radius` | z_B, B, Te_interp primitives; Te_J, r_L, r_L_ratio formula columns |
| `ExB_Drift` | z, E, B_on_cells primitives; v_ExB formula column |
| `Scalar_Metrics` | All scalar performance metrics as cross-sheet INDEX/MATCH/SUMPRODUCT formulas — no hardcoded values |

> Formula columns are labelled `[formula]`; primitive inputs have no suffix. Julia `[ref]` columns are provided for cross-check — values should match.

---

## Console Output

```
── DERIVED PERFORMANCE METRICS — MyHET ──────────────────────
  Δφ (acceleration potential) :  XXX.X  V
  v_exit (S_iz-weighted)      :   XX.XX  km/s
  Isp    (S_iz-weighted)      :  XXXX   s
  Thrust estimate             :    X.XX  mN
  η_m  (propellant util.)     :   XX.X  %
  η_v  (beam voltage frac.)   :   XX.X  %
  η_beam (η_m × η_v)          :   XX.X  %
  Power (beam estimate)       :  XXX.X  W
  E_max                       :  XXX.X  kV/m  at z = XX.X mm
  S_iz peak                   : X.XXe+XX  m⁻³s⁻¹  at z = XX.X mm
  Ionisation zone FWHM        :    X.X  mm
  Peak Hall parameter Ω_e     :   XXX   (at B_peak)
─────────────────────────────────────────────────────────────

Saved: MyHET_derived.png
Excel saved → MyHET_derived_data.xlsx
  Sheets: Constants | Electric_Field | Ion_Velocity | Ionisation_Source |
          Hall_Parameter | Larmor_Radius | ExB_Drift | Scalar_Metrics
```

---

## Dependencies

Inherited from `My_HET.jl` and `2D axisymmetric cross section.jl`.

| Package | Purpose |
|---|---|
| Plots.jl (GR backend) | All 6 panels |
| Plots.PlotMeasures | Margin control |
| XLSX.jl | 8-sheet Excel export with live formulas |
| Printf | Formatted console output |
| Statistics | Mean for beam current density estimate |

| Variable | Source |
|---|---|
| `viz_data` | `2D axisymmetric cross section.jl` |
| `z_peak` | `My_HET.jl` |

---

## References

- Goebel, D. M. & Katz, I. (2008). *Fundamentals of Electric Propulsion: Ion and Hall Thrusters*. JPL Space Science and Technology Series. Wiley. *(k_iz piecewise fit — Appendix E, Table E-1)*
- Mikellides, I., Katz, I. & Mandell, M. (2001). *A 1-D Model of the Hall-Effect Thruster with an Exhaust Region*. AIAA-2001-3505. *(original k_iz cross-section data source)*
