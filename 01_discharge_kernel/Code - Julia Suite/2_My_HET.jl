using HallThruster
using Printf
using Plots
using Plots.PlotMeasures
using XLSX

gr()

# ═════════════════════════════════════════════════════════════════════════
# MyHET_7mN — 7 mN Hall Effect Thruster Simulation
# HallThruster.jl v0.23.2 | Julia 1.12.6
# ═════════════════════════════════════════════════════════════════════════
#                    CONFIGURATION  —  EDIT HERE                        
# All tuneable parameters live in this block.                           
# Geometry, B-field, nominal operating point, sweep ranges, and anomalous transport coefficients.                                 
# The ConfigV3 script provides the optimized outputs for the desired Thrust and Isp. Paste the values below
# ============================================================================

# ── Geometry ──────────────────────────────────────────────────────────────────
const w_new     = 0.0099   # channel width        [m]  — 10 mm
const ratio     = 0.3318    # w / R_mean           [—]  — keep 0.25–0.35
const L_w_ratio = 2.621     # channel length / w   [—]  — keep 2.0–3.0

# ── Magnetic field ────────────────────────────────────────────────────────────
# B_peak ∝ 1/R_mean; anchored at reference thruster (R=32.5mm → B=180 Gauss)
const B_ref      = 0.01844    # reference peak B field  [T]
const R_ref      = 0.0325   # reference mean radius   [m]
const z_peak_frac = 0.900     # B peak as fraction of L [—]  — ionization near exit
const sigma_frac  = 0.341     # Gaussian width / L      [—]

# ── Nominal operating point ───────────────────────────────────────────────────
const Vd_nom   = 218      # discharge voltage    [V]
const mdot_ref = 7.0211703e-07 # anode mass flow rate [kg/s]  — 0.68 mg/s

# ── Anomalous transport (TwoZoneBohm) ─────────────────────────────────────────
# α_inner: inner zone coefficient  |  α_outer: outer zone coefficient
# Lower values → less cross-field transport → higher efficiency (but may not converge)
const α_inner = 1/160
const α_outer = 1/16

# ── 2D sweep ranges ───────────────────────────────────────────────────────────
const Vd_sweep   = [200.0, 300.0, 400.0, 500.0, 600.0]            # [V]
const mdot_sweep = [0.5, 0.65, 1.0, 2.0, 3.0, 4.0, 5.0] .* 1e-6  # [kg/s]

# ── Simulation domain & grid ──────────────────────────────────────────────────
const domain_end = 0.05     # domain length [m]  — channel (2.5 cm) + plume
const n_cells    = 50       # EvenGrid cell count

# ── SimParams: normal vs lean (mdot < 1 mg/s) ────────────────────────────────
const dt_normal       = 2e-9    # timestep for ṁ ≥ 1 mg/s  [s]
const dt_lean         = 1e-9    # timestep for ṁ < 1 mg/s  [s]
const duration_normal = 2e-4   # run duration for ṁ ≥ 1 mg/s  [s]
const duration_lean   = 4e-4   # run duration for ṁ < 1 mg/s  [s]
const num_save_normal = 100     # saved snapshots for ṁ ≥ 1 mg/s
const num_save_lean   = 300     # saved snapshots for ṁ < 1 mg/s

# ── Sanity filter for 2D grid ─────────────────────────────────────────────────
const isp_lo = 600.0    # minimum physical Isp  [s]
const isp_hi = 3000.0   # maximum physical Isp  [s]
const T_lo   = 0.5e-3   # minimum physical thrust [N]

# ══════════════════════════════════════════════════════════════════════════════
#                       END OF CONFIGURATION
# NOTE: η_PPU (PPU efficiency) is a hardware property of the power processing
# unit, NOT a simulation parameter.  Typical range: 0.85–0.90.
# P_system = P_discharge / η_PPU  (computed in Section 11 as a range)
# ══════════════════════════════════════════════════════════════════════════════
# ── DESIGN NOTES ──────────────────────────────────────────────────────────────
#   Channel width       : w = R_outer - R_inner
#   Mean radius         : R_mean = w / ratio
#   Width ratio         : w / R_mean  (healthy: 0.25–0.35)
#   Length/width ratio  : L / w       (healthy: 2.0–3.0)
#
#   Scaling laws:
#     Isp  ∝  √Vd           → double Vd gives Isp × 1.41
#     T    ∝  ṁ · √Vd       → double ṁ  gives T  × 2
#     T/P  ∝  1/√Vd         → higher Vd lowers thrust-per-watt
#     B_peak ∝ 1/R_mean     → larger thruster needs lower B
#
#   Voltage scaling table (η_v ≈ 0.9, Xe):
#     150 V → Isp ≈  780 s
#     200 V → Isp ≈  900 s
#     250 V → Isp ≈ 1005 s  ← nominal
#     300 V → Isp ≈ 1100 s
#     400 V → Isp ≈ 1270 s
#     600 V → Isp ≈ 1560 s
# ══════════════════════════════════════════════════════════════════════════════

# ── 1. GEOMETRY ───────────────────────────────────────────────────────────────
R_mean_new = w_new / ratio

my_geometry = HallThruster.Geometry1D(inner_radius = R_mean_new - w_new/2,
    outer_radius   = R_mean_new + w_new/2, channel_length = w_new * L_w_ratio,)

# ── 2. MAGNETIC FIELD ─────────────────────────────────────────────────────────
L = my_geometry.channel_length
B_max = B_ref * (R_ref / R_mean_new)
z_peak = z_peak_frac * L 
sigma = sigma_frac * L 

z_vec = collect(range(0.0, domain_end, length=100))
B_vec = B_max .* exp.(.-(z_vec .- z_peak) .^ 2 ./ (2 * sigma^2))

my_Bfield_data = HallThruster.MagneticField(file="", z=z_vec, B=B_vec,)

# ── 3. THRUSTER ───────────────────────────────────────────────────────────────
my_thruster = HallThruster.Thruster(
    name="MyHET_7mN",
    geometry=my_geometry,
    magnetic_field=my_Bfield_data,
    shielded=false,
)

# ── 4. CONFIG ─────────────────────────────────────────────────────────────────
my_config = HallThruster.Config(
    thruster=my_thruster,
    discharge_voltage = Vd_nom,
    anode_mass_flow_rate = mdot_ref,
    domain=(0.0, domain_end),
    anom_model=HallThruster.TwoZoneBohm(α_inner, α_outer),
)

my_sim = HallThruster.SimParams(
    grid = HallThruster.EvenGrid(n_cells),
    num_save = num_save_lean,
    dt = dt_lean,
    duration = duration_lean,
)

# ── 5. RUN NOMINAL ────────────────────────────────────────────────────────────────────
@printf("Running nominal simulation (%.0f V, %.4f mg/s)...\n", Vd_nom, mdot_ref*1e6)
@time sol = HallThruster.run_simulation(my_config, my_sim)
println("Retcode : ", sol.retcode)

# ── 6. NOMINAL RESULTS ────────────────────────────────────────────────────────────────
sol_avg = HallThruster.time_average(sol)

T   = HallThruster.thrust(sol_avg)[]
Id  = HallThruster.discharge_current(sol_avg)[]
η   = HallThruster.anode_eff(sol_avg)[]
Isp = T / (mdot_ref * 9.81)
P   = Vd_nom * Id

@printf("\n── Nominal Results (%.0f V, %.4f mg/s) ──────────\n", Vd_nom, mdot_ref*1e6)
@printf("Thrust            : %.2f mN\n", T * 1e3)
@printf("Discharge current : %.2f A\n",  Id)
@printf("Power             : %.1f W\n",  P)
@printf("Isp               : %.0f s\n",  Isp)
@printf("Anode efficiency  : %.1f %%\n", η * 100)

# ── 7. 2D PERFORMANCE GRID (Vd x ṁ)  ───────────────────────────────────────────────────────────
println("\n2D Performance Grid (Vd × ṁ)")
@printf("%-8s  %-9s  %-10s  %-8s  %-9s  %-9s\n",
    "Vd(V)", "ṁ(mg/s)", "T(mN)", "Isp(s)", "P(W)", "η(%)")
println(repeat("─", 62))

grid_results = []

for Vd in Vd_sweep, mdot in mdot_sweep

    is_lean = mdot < 1.0e-6

    cfg = HallThruster.Config(
        thruster             = my_thruster,
        discharge_voltage    = Vd,
        anode_mass_flow_rate = mdot,
        domain               = (0.0, domain_end),
        anom_model           = HallThruster.TwoZoneBohm(α_inner, α_outer),
    )
    sim = HallThruster.SimParams(
        grid     = HallThruster.EvenGrid(n_cells),
        num_save = is_lean ? num_save_lean   : num_save_normal,
        dt       = is_lean ? dt_lean         : dt_normal,
        duration = is_lean ? duration_lean   : duration_normal,
    )

    s = HallThruster.run_simulation(cfg, sim)

    if s.retcode == :success
        avg    = HallThruster.time_average(s)
        T_sw   = HallThruster.thrust(avg)[]
        Id_sw  = HallThruster.discharge_current(avg)[]
        η_sw   = HallThruster.anode_eff(avg)[]
        Isp_sw = T_sw / (mdot * 9.81)
        P_sw   = Vd * Id_sw

        if T_sw < T_lo || Isp_sw < isp_lo || Isp_sw > isp_hi
            @printf("%-8.0f  %-9.3f  UNPHYSICAL (T=%.3f mN, Isp=%.0f s) — skipped\n",
                Vd, mdot*1e6, T_sw*1e3, Isp_sw)
        else
            @printf("%-8.0f  %-9.3f  %-10.2f  %-8.0f  %-9.1f  %-9.1f\n",
                Vd, mdot*1e6, T_sw*1e3, Isp_sw, P_sw, η_sw*100)
            push!(grid_results, (Vd=Vd, mdot=mdot*1e6, T_mN=T_sw*1e3,
                                  Isp=Isp_sw, P_W=P_sw, eta=η_sw*100))
        end
    else
        @printf("%-8.0f  %-9.3f  FAILED (%s)\n", Vd, mdot*1e6, s.retcode)
    end

end

# ── 9. PLASMA PROFILES (Nominal Point) ───────────────────────────────────────────────────────────
z_cells = sol_avg.grid .* 100
exit_x = L*100

p1 = Plots.plot(z_vec .* 100, B_vec .* 1e4,
    xlabel="z (cm)", ylabel="B (Gauss)",
    title="Magnetic Field Profile", label="B(z)",
    color=:blue, linewidth=2,
    titlefontsize=10, guidefontsize=9, margin=6mm, legend=:topright,
)
Plots.vline!([exit_x], color=:red, linestyle=:dash, label="channel exit")

Te = sol_avg[:Tev][1]
p2 = Plots.plot(z_cells, Te,
    xlabel="z (cm)", ylabel="Te (eV)",
    title="Electron Temperature", label="Tev",
    color=:orange, linewidth=2,
    titlefontsize=10, guidefontsize=9, margin=6mm, legend=:topright,
)
Plots.vline!([exit_x], color=:red, linestyle=:dash, label=false)

phi = sol_avg[:potential][1]
p3 = Plots.plot(z_cells, phi,
    xlabel="z (cm)", ylabel="Potential (V)",
    title="Plasma Potential", label="ϕ(z)",
    color=:green, linewidth=2,
    titlefontsize=10, guidefontsize=9, margin=6mm, legend=:topright,
)
Plots.vline!([exit_x], color=:red, linestyle=:dash, label=false)

ji = sol_avg[:ji][1]
p4 = Plots.plot(z_cells, ji,
    xlabel="z (cm)", ylabel="Ion current (A/m²)",
    title="Ion Current Density", label="ji(z)",
    color=:purple, linewidth=2,
    titlefontsize=10, guidefontsize=9, margin=6mm, legend=:bottomright,
)
Plots.vline!([exit_x], color=:red, linestyle=:dash, label=false)

fig_plasma = Plots.plot(p1, p2, p3, p4,
    layout=(2, 2), size=(1400, 900),
    plot_title="MyHET_7mN — Plasma Properties ($(Int(Vd_nom)) V, $(round(mdot_ref*1e6, digits=4)) mg/s)",
    plot_titlefontsize=13,
    left_margin=12mm, right_margin=10mm, top_margin=10mm, bottom_margin=10mm,
)

savefig(fig_plasma, "MyHET_7mN_plasma.png")
println("Saved → MyHET_7mN_plasma.png")
display(fig_plasma)

# ── 9. PERFORMANCE CHARACTERIZATION PLOTS ──────────────────────────────────────────────────────────────────
# 9.1  T/P vs Isp
p_tp = Plots.plot(
    xlabel="Isp (s)", ylabel="T/P (mN/W)", title="Thrust/Power vs Isp",
    legend=:topright, titlefontsize=10, guidefontsize=9, margin=6mm,
)
for Vd in Vd_sweep
    rows = sort(filter(r -> r.Vd == Vd, grid_results), by=r -> r.Isp)
    isempty(rows) && continue
    Plots.plot!(p_tp,
        [r.Isp for r in rows], [r.T_mN/r.P_W for r in rows],
        label="$(Int(Vd)) V", linewidth=2, marker=:circle, markersize=4,
    )
end

# 9.2  Thrust vs mass flow
p_T = Plots.plot(
    xlabel="Mass flow rate ṁ (mg/s)", ylabel="Thrust (mN)", title="Thrust vs Mass Flow Rate",
    legend=:topleft, titlefontsize=10, guidefontsize=9, margin=6mm,
)
for Vd in Vd_sweep
    rows = sort(filter(r -> r.Vd == Vd, grid_results), by=r -> r.mdot)
    isempty(rows) && continue
    Plots.plot!(p_T,
        [r.mdot for r in rows], [r.T_mN for r in rows],
        label="$(Int(Vd)) V", linewidth=2, marker=:circle, markersize=4,
    )
end

# 9.3  Isp vs Vd + √Vd theory
p_isp = Plots.plot(
    xlabel="Discharge voltage Vd (V)", ylabel="Isp (s)", title="Isp vs Discharge Voltage",
    legend=:topleft, titlefontsize=10, guidefontsize=9, margin=6mm,
)
for mdot in mdot_sweep
    rows = sort(
        filter(r -> isapprox(r.mdot, mdot*1e6, rtol=1e-3), grid_results),
        by=r -> r.Vd,
    )
    isempty(rows) && continue
    Plots.plot!(p_isp,
        [r.Vd for r in rows], [r.Isp for r in rows],
        label="$(round(mdot*1e6, digits=2)) mg/s",
        linewidth=1.5, marker=:square, markersize=4,
    )
end

# √Vd theory line — anchored at 300 V, 1.0 mg/s (best-converged reference)
ref_pts = filter(r -> r.Vd == 300.0 && isapprox(r.mdot, 1.0, rtol=0.01), grid_results)
isempty(ref_pts) && (ref_pts = filter(r -> r.Vd == 300.0 && isapprox(r.mdot, 0.65, rtol=0.01), grid_results))
Isp_anchor  = ref_pts[1].Isp
Vd_theory   = collect(range(200.0, 600.0, length=100))
Isp_theory  = Isp_anchor .* sqrt.(Vd_theory ./ 300.0)
Plots.plot!(p_isp, Vd_theory, Isp_theory,
    label="√Vd theory (anchored 300 V / 1 mg/s)",
    linewidth=2, linestyle=:dash, color=:black,
)

fig_perf = Plots.plot(p_tp, p_T, p_isp,
    layout=(1, 3), size=(1600, 500),
    plot_title="MyHET_7mN — Performance Characterisation",
    plot_titlefontsize=12,
    left_margin=12mm, right_margin=8mm, top_margin=10mm, bottom_margin=10mm,
)

savefig(fig_perf, "MyHET_7mN_performance.png")
println("Saved → MyHET_7mN_performance.png")
display(fig_perf)

# ── 10. EXPORT TO EXCEL ───────────────────────────────────────────────────────
# Nominal summary pulled from grid_results at (Vd_nom, mdot_ref) for consistency.
# Falls back to standalone sol if that point was filtered as unphysical.
nom_grid = filter(r -> r.Vd == Vd_nom && isapprox(r.mdot, mdot_ref*1e6, rtol=0.01),
                  grid_results)
if !isempty(nom_grid)
    nom        = nom_grid[1]
    T_nom      = nom.T_mN * 1e-3
    Id_nom     = nom.P_W / Vd_nom
    η_nom      = nom.eta / 100.0
    Isp_nom_xl = nom.Isp
    P_nom      = nom.P_W
else
    T_nom      = T
    Id_nom     = Id
    η_nom      = η
    Isp_nom_xl = Isp
    P_nom      = P
end

XLSX.openxlsx("MyHET_7mN_results.xlsx", mode="w") do xf

    # Sheet 1: Plasma Properties
    sheet1 = xf[1]
    XLSX.rename!(sheet1, "Plasma Properties")
    sheet1["A1"] = "z (cm)"
    sheet1["B1"] = "B field (Gauss)"
    sheet1["C1"] = "Electron Temp Te (eV)"
    sheet1["D1"] = "Plasma Potential (V)"
    sheet1["E1"] = "Ion Current Density (A/m²)"

    for i in eachindex(z_cells)
        sheet1[i+1, 1] = z_cells[i]
        sheet1[i+1, 2] = B_max * exp(-((z_cells[i]/100 - z_peak)^2) / (2*sigma^2)) * 1e4
        sheet1[i+1, 3] = Te[i]
        sheet1[i+1, 4] = phi[i]
        sheet1[i+1, 5] = ji[i]
    end

    # Sheet 2: Performance Summary
    XLSX.addsheet!(xf, "Performance Summary")
    sheet2 = xf[2]
    sheet2["A1"] = "Parameter"
    sheet2["B1"] = "Value"
    sheet2["C1"] = "Unit"

    summary_data = [
        "Thrust"            T_nom*1e3                       "mN"
        "Discharge Current" Id_nom                          "A"
        "Power"             P_nom                           "W"
        "Isp"               Isp_nom_xl                      "s"
        "Anode Efficiency"  η_nom*100                       "%"
        "Discharge Voltage" Vd_nom                          "V"
        "Mass Flow Rate"    mdot_ref*1e6                    "mg/s"
        "Peak B Field"      B_max*1e4                       "Gauss"
        "B Peak Location"   z_peak*100                      "cm"
        "Channel Length"    L*100                           "cm"
        "Inner Radius"      my_geometry.inner_radius*1e3    "mm"
        "Outer Radius"      my_geometry.outer_radius*1e3    "mm"
        "Mean Radius"       R_mean_new*1e3                  "mm"
        "Channel Width"     w_new*1e3                       "mm"
        "w/R_mean"          ratio                           "—"
        "L/w"               L_w_ratio                       "—"
    ]

    for i in 1:size(summary_data, 1)
        sheet2[i+1, 1] = summary_data[i, 1]
        sheet2[i+1, 2] = summary_data[i, 2]
        sheet2[i+1, 3] = summary_data[i, 3]
    end

    # Sheet 3: 2D Performance Grid
    XLSX.addsheet!(xf, "2D Performance Grid")
    sheet3 = xf[3]
    sheet3["A1"] = "Discharge Voltage (V)"
    sheet3["B1"] = "Mass Flow (mg/s)"
    sheet3["C1"] = "Thrust (mN)"
    sheet3["D1"] = "Isp (s)"
    sheet3["E1"] = "Power (W)"
    sheet3["F1"] = "Anode Efficiency (%)"
    sheet3["G1"] = "T/P (mN/W)"

    for (i, r) in enumerate(grid_results)
        sheet3[i+1, 1] = r.Vd
        sheet3[i+1, 2] = r.mdot
        sheet3[i+1, 3] = r.T_mN
        sheet3[i+1, 4] = r.Isp
        sheet3[i+1, 5] = r.P_W
        sheet3[i+1, 6] = r.eta
        sheet3[i+1, 7] = r.T_mN / r.P_W
    end

    # Sheet 4: Magnetic Field Profile
    XLSX.addsheet!(xf, "Magnetic Field")
    sheet4 = xf[4]
    sheet4["A1"] = "z (m)"
    sheet4["B1"] = "z (cm)"
    sheet4["C1"] = "B field (Tesla)"
    sheet4["D1"] = "B field (Gauss)"

    for i in eachindex(z_vec)
        sheet4[i+1, 1] = z_vec[i]
        sheet4[i+1, 2] = z_vec[i] * 100
        sheet4[i+1, 3] = B_vec[i]
        sheet4[i+1, 4] = B_vec[i] * 1e4
    end

end

println("Excel saved → MyHET_7mN_results.xlsx")


# ── 11. NOMINAL SUMMARY PRINTOUT ─────────────────────────────────────────────
const g₀ = 9.80665

# Derived quantities for summary
Id_sum   = HallThruster.discharge_current(sol_avg)[]   # from solution, not P/V
η_c_sum  = HallThruster.current_eff(sol_avg)[]         # I_ion / I_d  (physical: 0–1)
TP_sum   = (T_nom * 1e3) / P_nom                       # mN/W
P_sys_lo = P_nom / 0.90                                # best-case PPU (η_PPU=90%)
P_sys_hi = P_nom / 0.85                                # worst-case PPU (η_PPU=85%)
w_ch     = my_geometry.outer_radius - my_geometry.inner_radius  # channel width
w_R_sum  = w_ch / R_mean_new                           # width ratio check
L_w_sum  = L / w_ch                                    # L/w ratio check

println("\n" * repeat("═", 56))
println("  MyHET_7mN — Nominal Point Summary")
println(repeat("═", 56))
println("  ── Operating Point ──────────────────────────────")
@printf("  Discharge voltage  :  %.0f V\n",          Vd_nom)
@printf("  Mass flow rate     :  %.4f mg/s\n",        mdot_ref*1e6)
println("  ── Performance ──────────────────────────────────")
@printf("  Thrust             :  %.2f mN\n",          T_nom*1e3)
@printf("  Isp                :  %.0f s\n",           Isp_nom_xl)
@printf("  Discharge power    :  %.1f W\n",           P_nom)
@printf("  System power (PPU) :  %.1f – %.1f W  (η_PPU=85–90%%)\n", P_sys_lo, P_sys_hi)
println("  ── Auxiliary power budget (steady-state, additional to PPU) ─")
println("     Cathode keeper/heater  :  2 – 5 W")
println("     Electromagnet coils    :  3 – 8 W")
println("     Flow controller/valves :  1 – 3 W")
@printf("     Total system (incl. aux):  %.0f – %.0f W\n", P_sys_lo + 6.0, P_sys_hi + 16.0)
@printf("  T/P                :  %.3f mN/W\n",        TP_sum)
@printf("  Anode efficiency   :  %.1f %%\n",          η_nom*100)
@printf("  Discharge current  :  %.3f A\n",           Id_sum)
@printf("  Current utilisation:  %.3f  (η_c = I_ion/I_d, healthy: 0.6–0.9)\n", η_c_sum)
println("  ── Geometry ─────────────────────────────────────")
@printf("  Inner radius       :  %.3f mm\n",          my_geometry.inner_radius*1e3)
@printf("  Outer radius       :  %.3f mm\n",          my_geometry.outer_radius*1e3)
@printf("  Mean radius        :  %.3f mm\n",          R_mean_new*1e3)
@printf("  Channel width      :  %.1f mm\n",          w_ch*1e3)
@printf("  Channel length     :  %.1f mm\n",          L*1e3)
@printf("  w/R_mean           :  %.3f  (healthy: 0.25–0.35)\n", w_R_sum)
@printf("  L/w                :  %.2f   (healthy: 2.0–3.0)\n",  L_w_sum)
println("  ── Magnetic Field ───────────────────────────────")
@printf("  B_peak             :  %.1f Gauss\n",       B_max*1e4)
@printf("  B peak location    :  %.1f mm from inlet\n", z_peak*1e3)
println(repeat("═", 56))
println("  Outputs:")
println("    MyHET_7mN_plasma.png      — plasma properties")
println("    MyHET_7mN_performance.png — T/P·Isp, T·ṁ, Isp·Vd")
println("    MyHET_7mN_results.xlsx    — all data, 4 sheets")
println(repeat("═", 56))
