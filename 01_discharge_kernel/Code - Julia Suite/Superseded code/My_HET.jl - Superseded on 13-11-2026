using HallThruster
using Printf
using Plots
using Plots.PlotMeasures
using XLSX

gr()

# ══════════════════════════════════════════════════════════════════════════════
# MyHET_7mN — 7 mN Hall Effect Thruster Simulation
# HallThruster.jl v0.23.2 | Julia 1.12.6
# ══════════════════════════════════════════════════════════════════════════════
#
# ── DESIGN SCALING RELATIONS ─────────────────────────────────────────────────
#
# PRELIMINARY DESIGN GEOMETRY VALUES
#
# my_geometry = HallThruster.Geometry1D(
#     inner_radius   = 0.0275,   # 27.5 mm |
#     outer_radius   = 0.0375,   # 37.5 mm | 
#     channel_length = 0.025     # 25 mm   |
# )
#
# GEOMETRY:
#   Channel width       : w = R_outer - R_inner = 37.5 - 27.5 = 10 mm
#   Mean radius         : R_mean = (R_outer + R_inner) / 2 = 32.5 mm
#   Width ratio         : w / R_mean = 10 / 32.5 ≈ 0.31   (healthy: 0.25–0.35)
#   Length/width ratio  : L / w = 25 / 10 = 2.5            (healthy: 2.0–3.0)
#   Channel area        : A = π(R_outer² - R_inner²) ≈ 20.4 cm²
#
#   Rules:
#     → Increasing L/w improves ionization but raises wall losses
#     → w/R_mean > 0.4 risks non-uniform discharge across channel width
#     → Scaling R_mean up 2× requires ṁ up ~4× for same thrust density
#
# DISCHARGE VOLTAGE:
#   Exhaust velocity    : v_ex = sqrt(2·e·η_v·V_d / m_Xe)   [m/s]
#   Isp                 : Isp  = v_ex / g0  ∝  sqrt(V_d)
#   Thrust              : T    = ṁ · v_ex   ∝  ṁ · sqrt(V_d)
#   Thrust/Power ratio  : T/P  ∝  1/sqrt(V_d)   ← decreases at high voltage
#
#   Voltage scaling table (Voltage utilization efficiency (η_v) ≈ 0.9, Xe propellant):
#     V_d = 150 V  →  Isp ≈  780 s
#     V_d = 200 V  →  Isp ≈  900 s
#     V_d = 300 V  →  Isp ≈ 1100 s   ← current operating point
#     V_d = 400 V  →  Isp ≈ 1270 s
#     V_d = 600 V  →  Isp ≈ 1560 s
#
#   Rules:
#     → Double V_d  : Isp × 1.41,  Thrust × 1.41,  Power × 2
#     → Double ṁ    : Thrust × 2,  Isp unchanged,   Power × 2
#
# MAGNETIC FIELD:
#   Hall parameter      : Ω = e·B / (m_e · ν_e)  >>  1  (electrons trapped)
#   Ion parameter       : Ω_i = e·B / (m_Xe · ν_i)  <<  1  (ions free)
#   Peak B scaling      : B_peak ∝ 1 / R_mean
#   At R_mean = 32.5 mm : B_peak = 180 Gauss  →  Ω ≈ 50  ✓
#
#   Rules:
#     → Larger thruster (bigger R_mean) → lower B_peak needed
#     → B_peak too low  : electrons not trapped, efficiency drops
#     → B_peak too high : excess wall losses, efficiency drops
#     → Peak location at 0.9×L places ionization zone near exit ✓
#
# QUICK PARAMETER CHANGE GUIDE:
#   Want more thrust?    → increase anode_mass_flow_rate (linear effect)
#   Want higher Isp?     → increase discharge_voltage (sqrt effect)
#   Want lower power?    → decrease both V_d and ṁ proportionally
#   Rescaling thruster?  → keep w/R_mean ≈ 0.31, L/w ≈ 2.5, B∝1/R_mean
# ══════════════════════════════════════════════════════════════════════════════

# ── 1. GEOMETRY ───────────────────────────────────────────────────────────────
# 15 mm width, scaled geometry
w_new = 0.010                    # 10 mm
ratio = 0.31                     # keep w/R_mean healthy (healthy: 0.25–0.35)
R_mean_new = w_new / ratio            # = 48.4 mm
L_w_ratio = 2.5                      # keep L/w = 2.5

my_geometry = HallThruster.Geometry1D(
    inner_radius=R_mean_new - w_new/2,   # 40.9 mm
    outer_radius=R_mean_new + w_new/2,   # 55.9 mm
    channel_length=w_new * L_w_ratio        # 37.5 mm
)

# Mass flow needs to scale with channel area to maintain thrust density
# A_new/A_old ≈ (55.9²-40.9²)/(37.5²-27.5²) ≈ 3.8×
# So bump ṁ up ~3–4× to maintain ionization quality

# ── 2. MAGNETIC FIELD ─────────────────────────────────────────────────────────
L = my_geometry.channel_length
# B_peak scales as 1/R_mean — reduce it
B_max = 0.018 * (0.0325 / R_mean_new)   # ≈ 0.012 T (120 Gauss)
#B_max  = 0.018     # 180 Gauss | scales as 1/R_mean — reduce if thruster grows
z_peak = 0.9 * L  # peak at 90% of channel length — ionization zone near exit
sigma = 0.3 * L  # Gaussian width — controls sharpness of B profile

z_vec = collect(range(0.0, 0.05, length=100))
B_vec = B_max .* exp.(.-(z_vec .- z_peak) .^ 2 ./ (2 * sigma^2))

my_Bfield_data = HallThruster.MagneticField(
    file="",
    z=z_vec,
    B=B_vec,
)

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
    discharge_voltage=300.0,       # V  | Isp ∝ sqrt(V_d), T ∝ sqrt(V_d)
    anode_mass_flow_rate=0.6500071e-6,# kg/s = 0.65 mg/s | T ∝ ṁ (linear)
    domain=(0.0, 0.05),
    anom_model=HallThruster.TwoZoneBohm(1/160, 1/16),
)

# ── 5. SIMPARAMS ──────────────────────────────────────────────────────────────
my_sim = HallThruster.SimParams(
    grid=HallThruster.EvenGrid(50),
    num_save=200,
    dt=2e-9,
    duration=2e-4,
)

# ── 6. RUN ────────────────────────────────────────────────────────────────────
println("Running simulation...")
@time sol = HallThruster.run_simulation(my_config, my_sim)
println("Retcode : ", sol.retcode)

# ── 7. RESULTS ────────────────────────────────────────────────────────────────
sol_avg = HallThruster.time_average(sol)

T = HallThruster.thrust(sol_avg)[]
Id = HallThruster.discharge_current(sol_avg)[]
η = HallThruster.anode_eff(sol_avg)[]
Isp = T / (0.6500071e-6 * 9.81)

@printf("\n── Results ───────────────────\n")
@printf("Thrust            : %.2f mN\n", T * 1e3)
@printf("Discharge current : %.2f A\n", Id)
@printf("Power             : %.1f W\n", 300.0 * Id)
@printf("Isp               : %.0f s\n", Isp)
@printf("Anode efficiency  : %.1f %%\n", η * 100)

# ── 8. THRUST SWEEP ───────────────────────────────────────────────────────────
println("\nThrust sweep — varying mass flow rate:")
@printf("%-12s %-12s %-12s %-10s %-12s\n",
    "mdot(mg/s)", "Thrust(mN)", "Isp(s)", "Power(W)", "η_anode(%)")

for mdot in [0.5, 0.6500071, 1.0, 2.0, 3.0, 4.0, 5.0] .* 1e-6
    cfg = HallThruster.Config(
        thruster=my_thruster,
        discharge_voltage=300.0,
        anode_mass_flow_rate=mdot,
        domain=(0.0, 0.05),
        anom_model=HallThruster.TwoZoneBohm(1/160, 1/16),
    )
    sim = HallThruster.SimParams(
        grid=HallThruster.EvenGrid(50),
        num_save=100,
        dt=2e-9,
        duration=2e-4,
    )
    s = HallThruster.run_simulation(cfg, sim)
    if s.retcode == :success
        avg = HallThruster.time_average(s)
        T_sw = HallThruster.thrust(avg)[]
        Id_sw = HallThruster.discharge_current(avg)[]
        η_sw = HallThruster.anode_eff(avg)[]
        Isp_sw = T_sw / (mdot * 9.81)
        @printf("%-12.2f %-12.2f %-12.0f %-10.1f %-12.1f\n",
            mdot*1e6, T_sw*1e3, Isp_sw, 300.0*Id_sw, η_sw*100)
    else
        @printf("%-12.2f  retcode: %s\n", mdot*1e6, s.retcode)
    end
end

# ── 9. PLOTS ──────────────────────────────────────────────────────────────────
z_cells = sol_avg.grid .* 100
exit_x = L * 100

p1 = Plots.plot(z_vec .* 100, B_vec .* 1e4,
    xlabel="z (cm)", ylabel="B (Gauss)",
    title="Magnetic Field Profile", label="B(z)",
    color=:blue, linewidth=2,
    titlefontsize=10, guidefontsize=9,
    margin=6mm, legend=:topright,
)
Plots.vline!([exit_x], color=:red, linestyle=:dash, label=false)

Te = sol_avg[:Tev][1]
p2 = Plots.plot(z_cells, Te,
    xlabel="z (cm)", ylabel="Te (eV)",
    title="Electron Temperature", label="Tev",
    color=:orange, linewidth=2,
    titlefontsize=10, guidefontsize=9,
    margin=6mm, legend=:topright,
)
Plots.vline!([exit_x], color=:red, linestyle=:dash, label=false)

phi = sol_avg[:potential][1]
p3 = Plots.plot(z_cells, phi,
    xlabel="z (cm)", ylabel="Potential (V)",
    title="Plasma Potential", label="ϕ(z)",
    color=:green, linewidth=2,
    titlefontsize=10, guidefontsize=9,
    margin=6mm, legend=:topright,
)
Plots.vline!([exit_x], color=:red, linestyle=:dash, label=false)

ji = sol_avg[:ji][1]
p4 = Plots.plot(z_cells, ji,
    xlabel="z (cm)", ylabel="Ion current (A/m²)",
    title="Ion Current Density", label="ji(z)",
    color=:purple, linewidth=2,
    titlefontsize=10, guidefontsize=9,
    margin=6mm, legend=:bottomright,
)
Plots.vline!([exit_x], color=:red, linestyle=:dash, label=false)

fig = Plots.plot(p1, p2, p3, p4,
    layout=(2, 2),
    size=(1400, 900),
    plot_title="MyHET_7mN — 7 mN Hall Thruster Simulation",
    plot_titlefontsize=13,
    left_margin=12mm,
    right_margin=10mm,
    top_margin=10mm,
    bottom_margin=10mm,
)

savefig(fig, "MyHET_7mN_results.png")
println("Saved → MyHET_7mN_results.png")
display(fig)

# ── 10. EXPORT TO EXCEL ───────────────────────────────────────────────────────
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

    data = [
        "Thrust" T*1e3 "mN"
        "Discharge Current" Id "A"
        "Power" 300.0*Id "W"
        "Isp" Isp "s"
        "Anode Efficiency" η*100 "%"
        "Discharge Voltage" 300.0 "V"
        "Mass Flow Rate" 0.65e-6*1e6 "mg/s"
        "Peak B Field" B_max*1e4 "Gauss"
        "B Peak Location" z_peak*100 "cm"
        "Channel Length" L*100 "cm"
        "Inner Radius" 27.5 "mm"
        "Outer Radius" 37.5 "mm"
    ]

    for i in 1:size(data, 1)
        sheet2[i+1, 1] = data[i, 1]
        sheet2[i+1, 2] = data[i, 2]
        sheet2[i+1, 3] = data[i, 3]
    end

    # Sheet 3: Thrust Sweep
    XLSX.addsheet!(xf, "Thrust Sweep")
    sheet3 = xf[3]
    sheet3["A1"] = "Mass Flow (mg/s)"
    sheet3["B1"] = "Thrust (mN)"
    sheet3["C1"] = "Isp (s)"
    sheet3["D1"] = "Power (W)"
    sheet3["E1"] = "Anode Efficiency (%)"

    sweep_mdots = [0.5, 0.65, 1.0, 2.0, 3.0, 4.0, 5.0] .* 1e-6
    row = 2
    for mdot in sweep_mdots
        cfg = HallThruster.Config(
            thruster=my_thruster,
            discharge_voltage=300.0,
            anode_mass_flow_rate=mdot,
            domain=(0.0, 0.05),
            anom_model=HallThruster.TwoZoneBohm(1/160, 1/16),
        )
        sim = HallThruster.SimParams(
            grid=HallThruster.EvenGrid(50), num_save=100, dt=2e-9, duration=2e-4,
        )
        s = HallThruster.run_simulation(cfg, sim)
        if s.retcode == :success
            avg = HallThruster.time_average(s)
            T_sw = HallThruster.thrust(avg)[]
            Id_sw = HallThruster.discharge_current(avg)[]
            η_sw = HallThruster.anode_eff(avg)[]
            Isp_sw = T_sw / (mdot * 9.81)
            sheet3[row, 1] = mdot * 1e6
            sheet3[row, 2] = T_sw * 1e3
            sheet3[row, 3] = Isp_sw
            sheet3[row, 4] = 300.0 * Id_sw
            sheet3[row, 5] = η_sw * 100
            row += 1
        end
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
