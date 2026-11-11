# My_HET.jl — Hall Effect Thruster Simulation

A 1D Hall Effect Thruster simulation using [HallThruster.jl](https://github.com/UM-PEPL/HallThruster.jl), targeting a **7 mN thrust-class** xenon thruster at 300 V discharge.

Covers geometry definition, Gaussian magnetic field profiling, a baseline simulation run, a mass-flow parametric sweep, plasma profile plots, and Excel export.

**HallThruster.jl v0.23.2 · Julia 1.12.6**

---

## Table of Contents

- [Thruster Specification](#thruster-specification)
- [Design Scaling Relations](#design-scaling-relations)
- [Usage](#usage)
- [Code Walkthrough](#code-walkthrough)
- [Outputs](#outputs)

---

## Thruster Specification

| Parameter | Value |
|---|---|
| Inner radius | 27.5 mm |
| Outer radius | 37.5 mm |
| Channel length | 25 mm |
| Channel width `w` | 10 mm |
| Mean radius `R_mean` | 32.5 mm |
| Width ratio `w/R_mean` | 0.31 |
| Length/width ratio `L/w` | 2.5 |
| Discharge voltage | 300 V |
| Nominal mass flow rate | 0.65 mg/s |
| Peak magnetic field | ≈ 120 Gauss (scaled) |
| B-peak location | 0.9 × L (near exit) |
| Anomalous transport model | TwoZoneBohm (1/160, 1/16) |
| Grid cells | 50 (even) |
| Timestep | 2 ns |
| Simulation duration | 0.2 ms |

---

## Design Scaling Relations

### Geometry

| Derived quantity | Formula | Value | Healthy range |
|---|---|---|---|
| Width ratio | `w / R_mean` | 0.31 | 0.25 – 0.35 |
| Length/width ratio | `L / w` | 2.5 | 2.0 – 3.0 |
| Channel area | `π(R_outer² − R_inner²)` | ≈ 20.4 cm² | — |

- Increasing `L/w` improves ionisation but raises wall losses
- `w/R_mean > 0.4` risks non-uniform discharge across the channel width
- Scaling `R_mean` up 2× requires `ṁ` up ~4× for the same thrust density

### Discharge Voltage

| Quantity | Relation |
|---|---|
| Exhaust velocity | `v_ex = √(2·e·η_v·V_d / m_Xe)` |
| Specific impulse | `Isp ∝ √V_d` |
| Thrust | `T ∝ ṁ · √V_d` |
| Thrust/power ratio | `T/P ∝ 1/√V_d` — decreases at high voltage |

Voltage scaling table (η_v ≈ 0.9, Xe propellant):

| V_d | Approx. Isp |
|---|---|
| 150 V | ≈ 780 s |
| 200 V | ≈ 900 s |
| **300 V** | **≈ 1100 s** ← nominal |
| 400 V | ≈ 1270 s |
| 600 V | ≈ 1560 s |

### Magnetic Field

| Quantity | Note |
|---|---|
| Hall parameter | `Ω = e·B / (m_e·ν_e) >> 1` — electrons trapped |
| Ion parameter | `Ω_i = e·B / (m_Xe·ν_i) << 1` — ions unmagnetised |
| Peak B scaling | `B_peak ∝ 1 / R_mean` |

At `R_mean = 32.5 mm`: `B_peak = 180 Gauss` → `Ω ≈ 50` ✓

### Quick Parameter Change Guide

| Goal | Action |
|---|---|
| More thrust | Increase `anode_mass_flow_rate` (linear) |
| Higher Isp | Increase `discharge_voltage` (sqrt effect) |
| Lower power | Decrease both `V_d` and `ṁ` proportionally |
| Rescale thruster | Keep `w/R_mean ≈ 0.31`, `L/w ≈ 2.5`, `B ∝ 1/R_mean` |

---

## Usage

### Dependencies

```julia
using HallThruster
using Printf
using Plots
using Plots.PlotMeasures
using XLSX
```

Install packages via:

```julia
using Pkg
Pkg.add(["HallThruster", "Plots", "XLSX"])
```

### Run

```bash
julia MyHET_7mN.jl
```

Outputs:
- `MyHET_7mN_results.png` — 2×2 plasma profile plots
- `MyHET_7mN_results.xlsx` — four-sheet Excel workbook

---

## Code Walkthrough

### 1. Geometry

Scaled channel geometry maintaining healthy width and length ratios.

```julia
w_new      = 0.010                        # 10 mm channel width
ratio      = 0.31                         # w/R_mean (healthy: 0.25–0.35)
R_mean_new = w_new / ratio                # = 32.3 mm
L_w_ratio  = 2.5                          # L/w = 2.5

my_geometry = HallThruster.Geometry1D(
    inner_radius   = R_mean_new - w_new/2,  # 27.3 mm
    outer_radius   = R_mean_new + w_new/2,  # 37.3 mm
    channel_length = w_new * L_w_ratio      # 25 mm
)
```

> Channel area scales as `A_new/A_old ≈ 3.8×` relative to the preliminary design.
> Mass flow rate should be increased ~3–4× to maintain ionisation quality at the new size.

### 2. Magnetic Field

Gaussian profile with peak near channel exit; `B_peak` scaled as `1/R_mean`.

```julia
L      = my_geometry.channel_length
B_max  = 0.018 * (0.0325 / R_mean_new)  # ≈ 0.012 T (120 Gauss)
z_peak = 0.9 * L                         # ionisation zone near exit
sigma  = 0.3 * L                         # Gaussian width

z_vec = collect(range(0.0, 0.05, length = 100))
B_vec = B_max .* exp.(.-(z_vec .- z_peak).^2 ./ (2 * sigma^2))

my_Bfield_data = HallThruster.MagneticField(file = "", z = z_vec, B = B_vec)
```

### 3. Thruster

```julia
my_thruster = HallThruster.Thruster(
    name           = "MyHET_7mN",
    geometry       = my_geometry,
    magnetic_field = my_Bfield_data,
    shielded       = false,
)
```

### 4. Configuration

```julia
my_config = HallThruster.Config(
    thruster             = my_thruster,
    discharge_voltage    = 300.0,         # V   | Isp ∝ √V_d
    anode_mass_flow_rate = 0.6500071e-6,  # kg/s = 0.65 mg/s | T ∝ ṁ
    domain               = (0.0, 0.05),
    anom_model           = HallThruster.TwoZoneBohm(1/160, 1/16),
)
```

### 5. Simulation Parameters

```julia
my_sim = HallThruster.SimParams(
    grid     = HallThruster.EvenGrid(50),
    num_save = 200,
    dt       = 2e-9,   # s
    duration = 2e-4,   # s
)
```

### 6. Run

```julia
println("Running simulation...")
@time sol = HallThruster.run_simulation(my_config, my_sim)
println("Retcode : ", sol.retcode)
```

### 7. Results

```julia
sol_avg = HallThruster.time_average(sol)

T   = HallThruster.thrust(sol_avg)[]
Id  = HallThruster.discharge_current(sol_avg)[]
η   = HallThruster.anode_eff(sol_avg)[]
Isp = T / (0.6500071e-6 * 9.81)

@printf("Thrust            : %.2f mN\n", T * 1e3)
@printf("Discharge current : %.2f A\n",  Id)
@printf("Power             : %.1f W\n",  300.0 * Id)
@printf("Isp               : %.0f s\n",  Isp)
@printf("Anode efficiency  : %.1f %%\n", η * 100)
```

### 8. Thrust Sweep

Parametric sweep across seven mass flow rates at fixed 300 V.

```julia
for mdot in [0.5, 0.6500071, 1.0, 2.0, 3.0, 4.0, 5.0] .* 1e-6
    cfg = HallThruster.Config(
        thruster             = my_thruster,
        discharge_voltage    = 300.0,
        anode_mass_flow_rate = mdot,
        domain               = (0.0, 0.05),
        anom_model           = HallThruster.TwoZoneBohm(1/160, 1/16),
    )
    sim = HallThruster.SimParams(
        grid = HallThruster.EvenGrid(50), num_save = 100, dt = 2e-9, duration = 2e-4,
    )
    s = HallThruster.run_simulation(cfg, sim)
    if s.retcode == :success
        avg    = HallThruster.time_average(s)
        T_sw   = HallThruster.thrust(avg)[]
        Id_sw  = HallThruster.discharge_current(avg)[]
        η_sw   = HallThruster.anode_eff(avg)[]
        Isp_sw = T_sw / (mdot * 9.81)
        @printf("%-12.2f %-12.2f %-12.0f %-10.1f %-12.1f\n",
                mdot*1e6, T_sw*1e3, Isp_sw, 300.0*Id_sw, η_sw*100)
    end
end
```

### 9. Plots

2×2 figure: magnetic field · electron temperature · plasma potential · ion current density.
Red dashed line marks the channel exit plane. Saved to `MyHET_7mN_results.png`.

```julia
z_cells = sol_avg.grid .* 100
exit_x  = L * 100

p1 = Plots.plot(z_vec .* 100, B_vec .* 1e4,
    xlabel = "z (cm)", ylabel = "B (Gauss)", title = "Magnetic Field Profile",
    label = "B(z)", color = :blue, linewidth = 2, margin = 6mm)
Plots.vline!([exit_x], color = :red, linestyle = :dash, label = false)

p2 = Plots.plot(z_cells, sol_avg[:Tev][1],
    xlabel = "z (cm)", ylabel = "Te (eV)", title = "Electron Temperature",
    label = "Tev", color = :orange, linewidth = 2, margin = 6mm)
Plots.vline!([exit_x], color = :red, linestyle = :dash, label = false)

p3 = Plots.plot(z_cells, sol_avg[:potential][1],
    xlabel = "z (cm)", ylabel = "Potential (V)", title = "Plasma Potential",
    label = "ϕ(z)", color = :green, linewidth = 2, margin = 6mm)
Plots.vline!([exit_x], color = :red, linestyle = :dash, label = false)

p4 = Plots.plot(z_cells, sol_avg[:ji][1],
    xlabel = "z (cm)", ylabel = "Ion current (A/m²)", title = "Ion Current Density",
    label = "ji(z)", color = :purple, linewidth = 2, margin = 6mm)
Plots.vline!([exit_x], color = :red, linestyle = :dash, label = false)

fig = Plots.plot(p1, p2, p3, p4,
    layout = (2, 2), size = (1400, 900),
    plot_title = "MyHET_7mN — 7 mN Hall Thruster Simulation",
    left_margin = 12mm, right_margin = 10mm, top_margin = 10mm, bottom_margin = 10mm)

savefig(fig, "MyHET_7mN_results.png")
```

### 10. Excel Export

Writes `MyHET_7mN_results.xlsx` with four sheets:

| Sheet | Contents |
|---|---|
| Plasma Properties | `z`, `B`, `Te`, `φ`, `ji` axial profiles |
| Performance Summary | Thrust, Isp, power, efficiency, geometry, operating point |
| Thrust Sweep | Sweep results across 7 mass flow rates |
| Magnetic Field | Full `z`-`B` profile in both SI and CGS units |

```julia
XLSX.openxlsx("MyHET_7mN_results.xlsx", mode = "w") do xf

    # Sheet 1 — Plasma Properties
    sheet1 = xf[1]
    XLSX.rename!(sheet1, "Plasma Properties")
    sheet1["A1"] = "z (cm)";  sheet1["B1"] = "B field (Gauss)"
    sheet1["C1"] = "Electron Temp Te (eV)"
    sheet1["D1"] = "Plasma Potential (V)";  sheet1["E1"] = "Ion Current Density (A/m²)"
    for i in eachindex(z_cells)
        sheet1[i+1, 1] = z_cells[i]
        sheet1[i+1, 2] = B_max * exp(-((z_cells[i]/100 - z_peak)^2) / (2*sigma^2)) * 1e4
        sheet1[i+1, 3] = Te[i];  sheet1[i+1, 4] = phi[i];  sheet1[i+1, 5] = ji[i]
    end

    # Sheet 2 — Performance Summary
    XLSX.addsheet!(xf, "Performance Summary")
    sheet2 = xf[2]
    sheet2["A1"] = "Parameter";  sheet2["B1"] = "Value";  sheet2["C1"] = "Unit"
    data = [
        "Thrust"            T*1e3       "mN"
        "Discharge Current" Id          "A"
        "Power"             300.0*Id    "W"
        "Isp"               Isp         "s"
        "Anode Efficiency"  η*100       "%"
        "Discharge Voltage" 300.0       "V"
        "Mass Flow Rate"    0.65        "mg/s"
        "Peak B Field"      B_max*1e4   "Gauss"
        "B Peak Location"   z_peak*100  "cm"
        "Channel Length"    L*100       "cm"
        "Inner Radius"      27.5        "mm"
        "Outer Radius"      37.5        "mm"
    ]
    for i in 1:size(data, 1)
        sheet2[i+1, 1] = data[i, 1];  sheet2[i+1, 2] = data[i, 2];  sheet2[i+1, 3] = data[i, 3]
    end

    # Sheet 3 — Thrust Sweep
    XLSX.addsheet!(xf, "Thrust Sweep")
    sheet3 = xf[3]
    sheet3["A1"] = "Mass Flow (mg/s)";  sheet3["B1"] = "Thrust (mN)"
    sheet3["C1"] = "Isp (s)";  sheet3["D1"] = "Power (W)";  sheet3["E1"] = "Anode Efficiency (%)"
    row = 2
    for mdot in [0.5, 0.65, 1.0, 2.0, 3.0, 4.0, 5.0] .* 1e-6
        cfg = HallThruster.Config(thruster = my_thruster, discharge_voltage = 300.0,
            anode_mass_flow_rate = mdot, domain = (0.0, 0.05),
            anom_model = HallThruster.TwoZoneBohm(1/160, 1/16))
        sim = HallThruster.SimParams(grid = HallThruster.EvenGrid(50),
            num_save = 100, dt = 2e-9, duration = 2e-4)
        s = HallThruster.run_simulation(cfg, sim)
        if s.retcode == :success
            avg = HallThruster.time_average(s)
            T_sw = HallThruster.thrust(avg)[];  Id_sw = HallThruster.discharge_current(avg)[]
            η_sw = HallThruster.anode_eff(avg)[];  Isp_sw = T_sw / (mdot * 9.81)
            sheet3[row, 1] = mdot*1e6;  sheet3[row, 2] = T_sw*1e3
            sheet3[row, 3] = Isp_sw;    sheet3[row, 4] = 300.0*Id_sw
            sheet3[row, 5] = η_sw*100;  row += 1
        end
    end

    # Sheet 4 — Magnetic Field Profile
    XLSX.addsheet!(xf, "Magnetic Field")
    sheet4 = xf[4]
    sheet4["A1"] = "z (m)";  sheet4["B1"] = "z (cm)"
    sheet4["C1"] = "B (Tesla)";  sheet4["D1"] = "B (Gauss)"
    for i in eachindex(z_vec)
        sheet4[i+1, 1] = z_vec[i];  sheet4[i+1, 2] = z_vec[i]*100
        sheet4[i+1, 3] = B_vec[i];  sheet4[i+1, 4] = B_vec[i]*1e4
    end

end
println("Excel saved → MyHET_7mN_results.xlsx")
```

---

## Outputs

| File | Description |
|---|---|
| `MyHET_7mN_results.png` | 2×2 axial plasma profile plots |
| `MyHET_7mN_results.xlsx` | Plasma properties · Performance summary · Thrust sweep · Magnetic field |

---

## References

- [HallThruster.jl documentation](https://um-pepl.github.io/HallThruster.jl/dev/)
- Rapp, D. & Francis, W. E. (1962). *Charge exchange between gaseous ions and atoms.* J. Chem. Phys., 37(11), 2631–2645.
