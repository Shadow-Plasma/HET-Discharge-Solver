# 2D Axisymmetric Cross Section — Data Extraction

Extracts and bundles plasma profile data from the 1D discharge simulation into a clean `viz_data` named tuple for downstream 2D visualisation.

> **Prerequisite:** Run `My_HET.jl` first in the same REPL session. This script reads variables already in memory — do not run it in a fresh REPL.

---

## Usage

```julia
# Step 1 — Run the discharge simulation first
julia> include("My_HET.jl")

# Step 2 — Extract and bundle data
julia> include("2D axisymmetric cross section.jl")
```

---

## What This Script Does

### 1. Geometry

Reads channel geometry directly from the `my_geometry` object in memory.

```julia
r_inner = my_geometry.inner_radius    # m
r_outer = my_geometry.outer_radius    # m
L_ch    = my_geometry.channel_length  # m
z_domain = my_config.domain           # (z_start, z_end) tuple
```

**Console output:**
```
── Geometry ──────────────────────────────
  Inner radius   : XX.X mm
  Outer radius   : XX.X mm
  Channel width  : XX.X mm
  Channel length : XX.X mm
  Domain         : 0 → 50 mm
```

---

### 2. Magnetic Field

Reads the axial B-field profile from `z_vec` and `B_vec` (set in `My_HET.jl`).

```julia
z_B    = z_vec    # axial positions (m), 100 points
B_prof = B_vec    # B magnitude (T), 100 points
```

**Console output:**
```
── Magnetic Field ────────────────────────
  B_max         : XXX.X Gauss  (X.XXXX T)
  Peak location : z = XX.XX mm  (XX.X% of channel)
  z_vec length  : 100 points
```

---

### 3. Plasma Profile Data

Extracts all plasma quantities from the time-averaged solution `sol_avg`.

```julia
z_cells  = sol_avg.grid          # axial cell centres (m)
Te_prof  = sol_avg[:Tev][1]      # electron temperature (eV)
phi_prof = sol_avg[:potential][1] # plasma potential (V)
ji_prof  = sol_avg[:ji][1]       # ion current density (A/m²)
ni_prof  = sol_avg[:ni][1]       # ion number density (m⁻³)
nn_prof  = sol_avg[:nn][1]       # neutral number density (m⁻³)
```

| Quantity | Symbol | Units |
|---|---|---|
| Electron temperature | `Te` | eV |
| Plasma potential | `φ` | V |
| Ion current density | `ji` | A/m² |
| Ion number density | `ni` | m⁻³ |
| Neutral number density | `nn` | m⁻³ |

---

### 4. Operating Conditions

```julia
Vd   = my_config.discharge_voltage   # V
mdot = 0.6500071e-6                  # kg/s (fallback if not in Config)
```

> **Note:** In HallThruster.jl v0.23.2 the mass flow rate is not directly accessible as a `Config` field. The script attempts `my_config.propellants[1].flow_rate` first and falls back to the session value `0.6500071e-6 kg/s` if unavailable.

---

### 5. viz_data Bundle

All extracted quantities are packed into a single named tuple for the 2D visualiser.

```julia
viz_data = (
    # Geometry
    r_inner, r_outer, L_ch, z_min, z_max,

    # B-field
    z_B, B_prof,

    # Plasma profiles
    z_cells, Te, phi, ji, ni, nn,

    # Operating conditions
    Vd, mdot_mgs,
)
```

**Console output:**
```
── viz_data bundle ready ─────────────────
  Fields: (:r_inner, :r_outer, :L_ch, :z_min, :z_max, :z_B, :B_prof,
           :z_cells, :Te, :phi, :ji, :ni, :nn, :Vd, :mdot_mgs)

Data extraction complete. Ready for visualizer.
```

---

## Output

No files are written. The script outputs `viz_data` into the REPL session scope for immediate use by `Visualise thruster 2d.jl`.

---

## Dependencies

Inherited from `My_HET.jl` — no additional packages required.

| Variable | Source |
|---|---|
| `my_geometry` | `My_HET.jl` |
| `my_config` | `My_HET.jl` |
| `sol_avg` | `My_HET.jl` |
| `z_vec`, `B_vec` | `My_HET.jl` |
