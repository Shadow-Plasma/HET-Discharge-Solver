# Visualise Thruster 2D — 2D Axisymmetric Visualiser

Generates an 8-panel 2D axisymmetric cross-section visualisation of the Hall thruster discharge and plume, using plasma profile data bundled by `2D axisymmetric cross section.jl`.

> **Prerequisites:** Run `My_HET.jl` and `2D axisymmetric cross section.jl` first in the same REPL session. This script reads `viz_data` from memory — it will fail in a fresh REPL.

---

## Usage

```julia
# Run all three in order in the same REPL session
julia> include("My_HET.jl")
julia> include("2D axisymmetric cross section.jl")
julia> include("Visualise thruster 2d.jl")
```

**Output file:** `MyHET_2D_viz.png` — saved in the current working directory.

---

## What This Script Does

### 0. Unpack viz_data

Reads all geometry, B-field, plasma profiles, and operating conditions from the `viz_data` named tuple. Converts units to mm and practical CGS/SI:

| Quantity | Conversion |
|---|---|
| Positions | m → mm |
| B-field | T → Gauss |
| Ion density `ni` | m⁻³ → ×10¹⁷ m⁻³ |
| Neutral density `nn` | m⁻³ → ×10¹⁹ m⁻³ |

---

### Zone Boundaries (Auto-Derived)

Zone boundaries are computed automatically from the geometry and B-field peak location — they update whenever the thruster is rescaled.

| Zone | Extent | Colour |
|---|---|---|
| Neutral Fill | 0 → 60% L_ch | — |
| Ionisation | 60% L_ch → z_B_peak | Yellow shading |
| Acceleration | z_B_peak → L_ch (exit) | Blue shading |
| Plume | L_ch → domain end | — |

---

### Plume Divergence Angle

Two divergence estimates are computed and shown on the geometry panel:

**Option 2 — ji 50%-decay half-angle (cyan cone):**

```
α_ji = arctan( r_c / Δz_50 )
```

Where `Δz_50` is the axial distance from the exit plane to where `ji` first reaches 50% of its plume plateau value. Near-field observable from the actual current density profile.

**Option 3 — Voltage utilisation (orange cone):**

```
α_volt = arccos( sqrt( 1 − Te_exit / Vd ) )
```

Geometry-independent. Uses only `Te_exit` and `Vd` from simulation outputs. Represents the ideal axial beam divergence from energy conservation.

> **Note:** Experimental plume angles are typically 2–3× larger than the voltage-utilisation value due to CEX ions, near-field measurement distances, and charge-state mixing.

Two additional estimates are printed to the console for reference:

| Method | Expression | Note |
|---|---|---|
| Option 1 — Geometric | `arctan(r_c / z_plume)` | Domain-dependent, reference only |
| Option 4 — φ residual | `arccos(√(1 − φ_exit/Vd))` | Potential-based, console reference |

---

### Panels

| Panel | Quantity | Colour |
|---|---|---|
| A | HET cross-section geometry schematic with zone shading and plume divergence cones | — |
| B | Magnetic field `B(z)` | Cyan |
| C | Electron temperature `Te(z)` | Orange |
| D | Plasma potential `φ(z)` | Green |
| E | Ion current density `ji(z)` | Purple |
| F | Ion number density `ni(z)` | Light blue |
| G | Neutral number density `nn(z)` | Pink-red |
| H | Normalised overlay — all profiles on the same 0–1 axis | Multi-colour |

All panels (except H) carry a secondary right Y-axis showing channel radial geometry reference lines at `±r_i` and `±r_o`.

---

### Panel A — Geometry Schematic

The main schematic panel includes:

- Channel walls (boron nitride) at `±r_i` and `±r_o`
- Inner core (anode region) shading
- Zone shading (ionisation, acceleration, plume)
- Exit plane dashed line
- B-peak location marker
- Dimension annotations (`r_o`, `r_i`, channel width `w`, `L_ch`)
- Two plume divergence cones (volt-util orange, ji-decay cyan)

---

### Console Output

```
Zone boundaries (auto-derived from geometry):
  Neutral Fill : 0  →  XX.X mm
  Ionisation   : XX.X  →  XX.X mm
  Acceleration : XX.X  →  XX.X mm (exit)
  Plume        : XX.X  →  XX.X mm

Plume divergence estimates:
  Exit-plane Te          : X.XX eV  (X.X% of Vd)
  Exit-plane phi         : X.XX V   (X.X% of Vd)
  ji plateau (plume)     : X.XX A/m²
  z at 50% ji plateau   : XX.X mm  (Δz=X.X mm from exit)
  ─────────────────────────────────────────────────────
  [PLOT] Volt-util  α    : XX.X °   ← Option 3
  [PLOT] ji-decay   α    : XX.X °   ← Option 2
  [REF]  Geometric  α    : XX.X °   ← Option 1 (domain-dependent)
  [REF]  phi-resid  α    : XX.X °   ← phi method
```

---

## Output

| File | Description |
|---|---|
| `MyHET_2D_viz.png` | 8-panel 2D axisymmetric visualisation, 1600×1800 px, 150 dpi |

---

## Dependencies

Inherited from `My_HET.jl` and `2D axisymmetric cross section.jl`.

| Package | Purpose |
|---|---|
| Plots.jl (GR backend) | All panels |
| Plots.PlotMeasures | Margin control |

| Variable | Source |
|---|---|
| `viz_data` | `2D axisymmetric cross section.jl` |
| `z_peak`, `L` | `My_HET.jl` |
