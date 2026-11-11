# Derived_Quantities.jl
# ─────────────────────────────────────────────────────────────────────────────
# Computes all derived plasma quantities from viz_data (no new simulation).
#
# Prerequisites
#   • My_HET.jl has been run
#   • 2D_axisymmetric_cross_section.jl has been run → viz_data in scope
#
# Usage
#   julia> include("Derived_Quantities.jl")
#
# Output
#   MyHET_derived.png   — 6-panel derived-quantities figure
# ─────────────────────────────────────────────────────────────────────────────

using Plots
using Plots.PlotMeasures
using Printf
using Statistics
gr()

# ── 0. Physical constants ─────────────────────────────────────────────────────
const e_q = 1.60218e-19   # C  (electron charge)
const me = 9.10938e-31   # kg (electron mass)
const Mi = 2.18014e-25   # kg (Natural Xe standard atomic weight, 131.293 u;
#     note: Xe-131 exact isotope = 130.905 u → 2.1737e-25 kg)
const g0 = 9.80665       # m/s²
const eps0 = 8.85419e-12   # F/m
const kB = 1.38065e-23   # J/K

# ── 1. Unpack viz_data (raw SI, no scaling) ───────────────────────────────────
r_i = viz_data.r_inner          # m
r_o = viz_data.r_outer          # m
L_ch = viz_data.L_ch             # m
z_min = viz_data.z_min            # m
z_max = viz_data.z_max            # m

z_B_m = vec(viz_data.z_B)         # m
B_T = vec(viz_data.B_prof)      # T

z_m = vec(viz_data.z_cells)     # m
Te_eV = vec(viz_data.Te)          # eV
phi_V = vec(viz_data.phi)         # V
ji_Am2 = vec(viz_data.ji)          # A/m²
ni_m3 = vec(viz_data.ni)          # m⁻³
nn_m3 = vec(viz_data.nn)          # m⁻³

Vd = viz_data.Vd               # V
mdot_mgs = viz_data.mdot_mgs         # mg/s
mdot_kgs = mdot_mgs * 1e-6           # kg/s

# mm versions for plotting x-axis
z_mm = z_m .* 1e3
z_B_mm = z_B_m .* 1e3
L_ch_mm = L_ch * 1e3

# Channel cross-section area (annular)
A_ch = π * (r_o^2 - r_i^2)          # m²

# ── 2. Numerical derivative helper (central difference, edge-padded) ──────────
function deriv(f, x)
    n = length(f)
    d = similar(f)
    for i in 2:(n-1)
        d[i] = (f[i+1] - f[i-1]) / (x[i+1] - x[i-1])
    end
    d[1] = (f[2] - f[1]) / (x[2] - x[1])
    d[n] = (f[n] - f[n-1]) / (x[n] - x[n-1])
    return d
end

# ── 3. Derived profiles ───────────────────────────────────────────────────────

# --- 3a. Electric field E(z) = -dφ/dz  [V/m] --------------------------------
E_field = -deriv(phi_V, z_m)          # V/m (positive = pointing toward anode)

# --- 3b. Ion velocity v_ion(z) from energy conservation ----------------------
# v²/2 = e·(φ_max - φ(z)) / Mi  — ions accelerated from high φ to low φ
phi_max = maximum(phi_V)
dV = phi_max .- phi_V           # potential drop experienced by ion
dV = max.(dV, 0.0)             # clip negative (numerical noise near anode)
v_ion = sqrt.(2 .* e_q .* dV ./ Mi)   # m/s

# --- 3c. Ionisation rate coefficient k_iz(Te) for Xenon ----------------------
# Piecewise fit from Goebel & Katz (2008), "Fundamentals of Electric Propulsion,"
# Appendix E — Table E-1.  Original data source:
#   Mikellides, Katz & Mandell, "A 1-D Model of the Hall-Effect Thruster with an
#   Exhaust Region," AIAA-2001-3505, 37th JPC, Salt Lake City, July 2001.
# Cross-sections measured for Xe; averaged over Maxwellian EEDF.
# Valid range: Te ≈ 0.5–10 eV (Table E-1). Accuracy: ~2–5% vs tabulated values.
#
#  Te < 5 eV:  ⟨σᵢvₑ⟩ = 1e-20·(3.97 + 0.643·Te − 0.0368·Te²)·exp(−12.127/Te)·v̄ₑ
#  Te ≥ 5 eV:  ⟨σᵢvₑ⟩ = 1e-20·(−1.031e-4·Te² + 6.386·exp(−12.127/Te))·v̄ₑ
#  where v̄ₑ = √(8·e·Te / π·mₑ)  [Maxwellian mean speed]
function k_iz_goebel(Te::Float64)::Float64
    Te = max(Te, 0.5)                                 # floor at 0.5 eV (formula lower bound)
    v_bar = sqrt(8.0 * e_q * Te / (π * me))           # mean electron speed [m/s]
    if Te < 5.0
        sigma_bar = 1e-20 * (3.97 + 0.643*Te - 0.0368*Te^2) * exp(-12.127/Te)
    else
        sigma_bar = 1e-20 * (-1.031e-4 * Te^2 + 6.386 * exp(-12.127/Te))
    end
    return sigma_bar * v_bar   # m³/s
end
k_iz = k_iz_goebel.(Te_eV)    # m³/s — element-wise over profile

# --- 3d. Ionisation source term S_iz(z)  [m⁻³ s⁻¹] -------------------------
# ne ≈ ni by quasineutrality
S_iz = k_iz .* nn_m3 .* ni_m3        # m⁻³ s⁻¹

# --- 3e. Hall parameter Ω_e(z) -----------------------------------------------
# Ω_e = e·B / (me·ν_eff)
# ν_eff estimated as electron-neutral collision freq + anomalous
# ν_en  ≈ k_en · nn,  k_en ≈ 2.5e-13 m³/s for Xe at ~10 eV
# This is order-of-magnitude; exact value needs cross-section data
k_en = 2.5e-13                        # m³/s  (e-n momentum transfer, Xe ~10 eV)
nu_en = k_en .* nn_m3                  # s⁻¹  (electron-neutral collision freq)
omega_ce = e_q .* B_T ./ me            # s⁻¹  (cyclotron frequency at B_field grid)
# Interpolate ν_en onto B_field z-grid for Ω_e
nu_en_interp = [nu_en[argmin(abs.(z_m .- zb))] for zb in z_B_m]
Omega_e = omega_ce ./ (nu_en_interp .+ 1.0)  # +1 avoids div-by-zero in plume

# --- 3f. Electron Larmor radius r_L(z) [mm] ----------------------------------
# r_L = sqrt(2·me·Te_J) / (e·B)  —  thermal Larmor radius
Te_J_interp = [Te_eV[argmin(abs.(z_m .- zb))] .* e_q for zb in z_B_m]
r_L_m = sqrt.(2 .* me .* Te_J_interp) ./ (e_q .* B_T .+ 1e-10)  # m
r_L_mm = r_L_m .* 1e3                 # mm
channel_width_mm = (r_o - r_i) * 1e3  # mm

# --- 3g. E×B drift velocity v_ExB(z) [km/s] ---------------------------------
# v_ExB = E / B — interpolate B onto cell grid
B_on_cells = [B_T[argmin(abs.(z_B_m .- zc))] for zc in z_m]
v_ExB = abs.(E_field) ./ (B_on_cells .+ 1e-10)   # m/s
v_ExB_kms = v_ExB ./ 1e3                          # km/s

# ── 4. Scalar performance metrics ─────────────────────────────────────────────
println("\n", "─"^52)
println("  DERIVED PERFORMANCE METRICS — MyHET")
println("─"^52)

# ── Exit / anode conditions ───────────────────────────────────────────────────
phi_exit_val = phi_V[end]                        # V  at last cell (plume)
phi_anode = maximum(phi_V)                    # V  peak potential (near anode)
delta_phi = phi_anode - phi_exit_val          # V  total acceleration potential

# ── S_iz-weighted mean exit velocity (corrected Isp) ─────────────────────────
# Each ion born at position z_i sees only (φ(z_i) − φ_exit) of the potential drop.
# Weight each birth location by the local ionisation source S_iz(z).
dV_local = max.(phi_V .- phi_exit_val, 0.0)         # V  per-cell potential drop
v_local = sqrt.(2.0 .* e_q .* dV_local ./ Mi)     # m/s  per-cell exit speed

S_on_cells = k_iz .* nn_m3 .* ni_m3                   # m⁻³s⁻¹  ionisation source
weights = S_on_cells ./ max(sum(S_on_cells), 1e-30) # normalised birth distribution

v_exit_ms = sum(weights .* v_local)                   # m/s  S_iz-weighted mean
v_exit_kms = v_exit_ms / 1e3
Isp_calc = v_exit_ms / g0                            # s

# ── Ion beam current from ji profile ─────────────────────────────────────────
# Use mean of the top-10% ji values as beam current density (robust vs single-cell spike)
ji_sorted = sort(abs.(ji_Am2), rev=true)
ji_top10 = mean(ji_sorted[1:max(1, div(length(ji_sorted), 10))])
I_beam = ji_top10 * A_ch                             # A
mdot_ion = I_beam * Mi / e_q                           # kg/s

# ── Efficiency metrics ────────────────────────────────────────────────────────
Thrust_est = mdot_ion * v_exit_ms                        # N
eta_m = min(mdot_ion / mdot_kgs, 1.0)              # propellant utilisation (cap at 100%)
eta_v = min(delta_phi / Vd, 1.0)                   # beam voltage fraction  (cap at 100%)
eta_beam = eta_m * eta_v                               # beam efficiency (no divergence in 1D)
Power_in = Vd * I_beam                                 # W

# ── E-field metrics ───────────────────────────────────────────────────────────
E_max = maximum(abs.(E_field))
z_Emax_mm = z_mm[argmax(abs.(E_field))]

# ── Ionisation zone metrics ───────────────────────────────────────────────────
S_peak = maximum(S_on_cells)
z_Speak = z_mm[argmax(S_on_cells)]
iz_half = findall(S_on_cells .>= 0.5 * S_peak)
iz_width = isempty(iz_half) ? NaN : z_mm[last(iz_half)] - z_mm[first(iz_half)]

# ── Print ─────────────────────────────────────────────────────────────────────
@printf "  Δφ (acceleration potential) : %6.1f  V\n" delta_phi
@printf "  v_exit (S_iz-weighted)      : %6.2f  km/s\n" v_exit_kms
@printf "  Isp    (S_iz-weighted)      : %6.0f  s\n" Isp_calc
@printf "  Thrust estimate             : %6.2f  mN\n" Thrust_est * 1e3
@printf "  η_m  (propellant util.)     : %6.1f  %%\n" eta_m * 100
@printf "  η_v  (beam voltage frac.)   : %6.1f  %%\n" eta_v * 100
@printf "  η_beam (η_m × η_v)          : %6.1f  %%\n" eta_beam * 100
@printf "  Power (beam estimate)       : %6.1f  W\n" Power_in
@printf "  E_max                       : %6.1f  kV/m  at z = %.1f mm\n" E_max/1e3 z_Emax_mm
@printf "  S_iz peak                   : %6.2e  m⁻³s⁻¹  at z = %.1f mm\n" S_peak z_Speak
@printf "  Ionisation zone FWHM        : %6.1f  mm\n" iz_width
@printf "  Peak Hall parameter Ω_e     : %6.0f  (at B_peak)\n" maximum(Omega_e)
println("─"^52)

# ── 5. Colour palette (matches main visualiser) ───────────────────────────────
BG = RGB(0.08, 0.08, 0.12)
C_Te = RGB(1.00, 0.55, 0.10)
C_phi = RGB(0.20, 0.85, 0.40)
C_ji = RGB(0.60, 0.30, 1.00)
C_ni = RGB(0.10, 0.75, 1.00)
C_nn = RGB(1.00, 0.30, 0.45)
C_B = RGB(0.50, 0.90, 1.00)
C_E = RGB(1.00, 0.85, 0.20)
C_v = RGB(0.80, 0.60, 1.00)
C_S = RGB(0.30, 1.00, 0.60)
C_Om = RGB(1.00, 0.60, 0.80)
C_rL = RGB(0.60, 0.90, 1.00)

_plt = (
    background_color=BG,
    background_color_inside=BG,
    foreground_color=:white,
    foreground_color_axis=:white,
    foreground_color_text=:white,
    foreground_color_grid=RGBA(1, 1, 1, 0.08),
    gridalpha=0.15,
    gridlinewidth=0.5,
    tickfontcolor=:white,
    guidefontcolor=:white,
    titlefontcolor=:white,
    legendfontcolor=:white,
    legendbackground_color=RGBA(0.08, 0.08, 0.12, 0.7),
    framestyle=:box,
    tick_direction=:out,
    minorticks=true,
    fontfamily="DejaVu Sans",
)

# ── Zone boundaries — derived from My_HET.jl geometry (auto-updates) ─────────
z_neutral_end_mm = 0.60 * L_ch_mm        # 60% of channel length
z_ionise_end_mm = z_peak * 1e3          # B-field peak location from My_HET.jl
z_accel_end_mm = L_ch_mm               # channel exit plane

function add_zones!(p)
    vspan!(p, [z_neutral_end_mm, z_ionise_end_mm];
        alpha=0.07, color=:yellow, label=false)
    vspan!(p, [z_ionise_end_mm, z_accel_end_mm];
        alpha=0.10, color=:deepskyblue, label=false)
    vline!(p, [L_ch_mm]; color=:white, linestyle=:dash, linewidth=1.0, alpha=0.5, label=false)
end

# ── 6. Panel 1 — Electric field E(z) ─────────────────────────────────────────
println("Building Panel 1: Electric field...")
p1 = plot(z_mm, E_field ./ 1e3;
    color=C_E, linewidth=2.0,
    xlabel="z (mm)", ylabel="E  (kV/m)",
    title="Electric Field  E(z) = −dφ/dz", titlefontsize=10,
    label="E(z)", _plt...)
hline!(p1, [0.0]; color=:white, linestyle=:dot, linewidth=0.6, alpha=0.3, label=false)
add_zones!(p1)
scatter!(p1, [z_Emax_mm], [E_max/1e3];
    color=C_E, markersize=5, markerstrokecolor=:white,
    label="E_max = $(round(E_max/1e3, digits=1)) kV/m")
annotate!(p1, z_Emax_mm+1.5, E_max/1e3*0.90,
    text("$(round(E_max/1e3,digits=1)) kV/m\n@ z=$(round(z_Emax_mm,digits=1)) mm",
        7, :white, :left))

# ── 7. Panel 2 — Ion velocity v_ion(z) ───────────────────────────────────────
println("Building Panel 2: Ion velocity...")
p2 = plot(z_mm, v_ion ./ 1e3;
    color=C_v, linewidth=2.0,
    xlabel="z (mm)", ylabel="v_ion  (km/s)",
    title="Ion Velocity  v(z) from Energy Conservation", titlefontsize=10,
    label="v_ion", _plt...)
add_zones!(p2)
annotate!(p2, 35.0, maximum(v_ion)/1e3 * 0.90,
    text("v_exit = $(round(v_exit_kms, digits=1)) km/s\nIsp = $(round(Isp_calc, digits=0)) s",
        8, :white, :left))

# ── 8. Panel 3 — Ionisation source S_iz(z) ───────────────────────────────────
println("Building Panel 3: Ionisation source...")
p3 = plot(z_mm, S_iz ./ 1e22;
    color=C_S, linewidth=2.0,
    xlabel="z (mm)", ylabel="S_iz  (×10²² m⁻³s⁻¹)",
    title="Ionisation Source  S_iz = k_iz(Te)·nn·ni", titlefontsize=10,
    label="S_iz", _plt...)
add_zones!(p3)
scatter!(p3, [z_Speak], [S_peak/1e22];
    color=C_S, markersize=5, markerstrokecolor=:white, label="Peak")
annotate!(p3, z_Speak+1.5, S_peak/1e22*0.88,
    text("Peak @ z=$(round(z_Speak,digits=1)) mm\nFWHM ≈ $(round(iz_width,digits=1)) mm",
        7, :white, :left))

# ── 9. Panel 4 — Hall parameter Ω_e(z) ───────────────────────────────────────
println("Building Panel 4: Hall parameter...")
p4 = plot(z_B_mm, Omega_e;
    color=C_Om, linewidth=2.0,
    xlabel="z (mm)", ylabel="Ω_e  = ω_ce / ν_en",
    title="Electron Hall Parameter  Ω_e(z)", titlefontsize=10,
    label="Ω_e", _plt...)
hline!(p4, [1.0]; color=:white, linestyle=:dash, linewidth=0.8, alpha=0.5,
    label="Ω_e = 1  (magnetised threshold)")
add_zones!(p4)
annotate!(p4, 28.0, maximum(Omega_e)*0.80,
    text("Max Ω_e = $(round(maximum(Omega_e), digits=0))", 8, :white, :left))

# ── 10. Panel 5 — Larmor radius r_L(z) / channel width ───────────────────────
println("Building Panel 5: Larmor radius...")
r_L_ratio = r_L_mm ./ channel_width_mm
p5 = plot(z_B_mm, r_L_mm;
    color=C_rL, linewidth=2.0,
    xlabel="z (mm)", ylabel="r_L  (mm)",
    title="Electron Larmor Radius  r_L(z)", titlefontsize=10,
    label="r_L", _plt...)
hline!(p5, [channel_width_mm];
    color=:yellow, linestyle=:dash, linewidth=1.0, alpha=0.7,
    label="Channel width = $(round(channel_width_mm, digits=1)) mm")
add_zones!(p5)
annotate!(p5, 30.0, channel_width_mm*1.05,
    text("r_L << channel → magnetised", 7, :yellow, :left))

# ── 11. Panel 6 — E×B drift velocity v_ExB(z) ────────────────────────────────
println("Building Panel 6: ExB drift...")
p6 = plot(z_mm, v_ExB_kms;
    color=C_B, linewidth=2.0,
    xlabel="z (mm)", ylabel="v_ExB  (km/s)",
    title="E×B Drift Velocity  v_ExB = E/B", titlefontsize=10,
    label="v_ExB", _plt...)
add_zones!(p6)
scatter!(p6, [z_mm[argmax(v_ExB_kms)]], [maximum(v_ExB_kms)];
    color=C_B, markersize=5, markerstrokecolor=:white, label="Peak drift")
annotate!(p6, z_mm[argmax(v_ExB_kms)]+1.5, maximum(v_ExB_kms)*0.90,
    text("$(round(maximum(v_ExB_kms),digits=0)) km/s", 7, :white, :left))

# ── 12. Compose & save ────────────────────────────────────────────────────────
println("Composing derived-quantities figure...")

fig2 = plot(
    p1, p2,
    p3, p4,
    p5, p6;
    layout=(3, 2),
    size=(1600, 1350),
    background_color=BG,
    plot_title="MyHET  —  Derived Plasma Quantities\n" *
        "Vd=$(Int(round(Vd))) V  |  mdot=$(round(mdot_mgs,digits=4)) mg/s  |  TwoZoneBohm(1/160,1/16)",
    plot_titlefontsize=11,
    plot_titlefontcolor=:white,
    left_margin=12mm,
    right_margin=10mm,
    top_margin=10mm,
    bottom_margin=10mm,
    dpi=150,
)

outfile2 = "MyHET_derived.png"
savefig(fig2, outfile2)
println("\nSaved: $(abspath(outfile2))")
display(fig2)

# ── 13. Excel export (live formulas) ─────────────────────────────────────────
# Design principle
#   • Each data sheet holds only the PRIMITIVE simulation outputs (the values
#     Julia computed that Excel cannot rederive from other columns).
#   • Every column that CAN be derived from primitives is an Excel formula,
#     so changing a primitive cell ripples through the whole sheet instantly.
#   • Scalar_Metrics sheet uses cross-sheet formulas referencing the data
#     sheets — no hardcoded numbers at all.
#
# Physical constants embedded in formulas  (same values as Julia block above)
#   e_q  = 1.60218e-19  C
#   me   = 9.10938e-31  kg
#   Mi   = 2.18014e-25  kg  (Natural Xe, 131.293 u)
#   g0   = 9.80665      m/s²
#   k_en = 2.5e-13      m³/s
#   k_iz = Goebel & Katz (2008) App. E piecewise fit — see Sheet 4 column E formula
# ─────────────────────────────────────────────────────────────────────────────

println("\nExporting data to Excel (with live formulas)...")
using XLSX

xlsx_file = "MyHET_derived_data.xlsx"

# How many data rows in each sheet (data starts at row 2)
N_cells = length(z_mm)   # 52  — plasma cell grid
N_B = length(z_B_mm) # 100 — B-field grid

XLSX.openxlsx(xlsx_file, mode="w") do xf

    # ════════════════════════════════════════════════════════════════════════
    # SHEET 1 — Constants  (named reference block used by formula sheets)
    # ════════════════════════════════════════════════════════════════════════
    sh1 = XLSX.addsheet!(xf, "Constants")
    sh1["A1"] = "Name"
    sh1["B1"] = "Value"
    sh1["C1"] = "Unit"
    sh1["D1"] = "Description"
    sh1["A2"] = "e_q"
    sh1["B2"] = 1.60218e-19
    sh1["C2"] = "C"
    sh1["D2"] = "Elementary charge"
    sh1["A3"] = "me"
    sh1["B3"] = 9.10938e-31
    sh1["C3"] = "kg"
    sh1["D3"] = "Electron mass"
    sh1["A4"] = "Mi"
    sh1["B4"] = 2.18014e-25
    sh1["C4"] = "kg"
    sh1["D4"] = "Natural Xe standard atomic weight (131.293 u)"
    sh1["A5"] = "g0"
    sh1["B5"] = 9.80665
    sh1["C5"] = "m/s²"
    sh1["D5"] = "Standard gravity"
    sh1["A6"] = "k_en"
    sh1["B6"] = 2.5e-13
    sh1["C6"] = "m³/s"
    sh1["D6"] = "e-n momentum transfer coeff (Xe ~10 eV)"
    sh1["A7"] = "kiz_Te_thr"
    sh1["B7"] = 5.0
    sh1["C7"] = "eV"
    sh1["D7"] = "Goebel k_iz piecewise threshold Te"
    sh1["A8"] = "kiz_E_thr"
    sh1["B8"] = 12.127
    sh1["C8"] = "eV"
    sh1["D8"] = "Goebel k_iz ionisation threshold (Xe IP in fit)"
    sh1["A9"] = "Vd"
    sh1["B9"] = Vd
    sh1["C9"] = "V"
    sh1["D9"] = "Discharge voltage"
    sh1["A10"] = "mdot_kgs"
    sh1["B10"] = mdot_kgs
    sh1["C10"] = "kg/s"
    sh1["D10"] = "Propellant mass flow rate"
    sh1["A11"] = "A_ch_m2"
    sh1["B11"] = A_ch
    sh1["C11"] = "m²"
    sh1["D11"] = "Channel cross-section area"
    sh1["A12"] = "ch_w_mm"
    sh1["B12"] = channel_width_mm
    sh1["C12"] = "mm"
    sh1["D12"] = "Channel width r_o - r_i"

    # ════════════════════════════════════════════════════════════════════════
    # SHEET 2 — Electric_Field
    # Primitives  : z_mm (A), phi_V (B)
    # Formula cols: E_kVm (C) = -(dφ/dz)/1000  via central difference in Excel
    #               [Note: central diff needs neighbours; write Julia value in D
    #                as a cross-check; E col is the formula version]
    # ════════════════════════════════════════════════════════════════════════
    sh2 = XLSX.addsheet!(xf, "Electric_Field")
    # Headers
    sh2["A1"] = "z_mm"
    sh2["B1"] = "phi_V"
    sh2["C1"] = "E_kVm [formula]"
    sh2["D1"] = "E_kVm [Julia ref]"
    sh2["E1"] = "Formula used"

    # Row 2 note
    sh2["E2"] = "C = -(phi[i+1]-phi[i-1])/(z[i+1]-z[i-1]) /1000  (central diff; edge rows use forward/backward)"

    for i in 1:N_cells
        r = i + 1   # Excel row
        sh2["A$r"] = z_mm[i]
        sh2["B$r"] = phi_V[i]
        # E = -dφ/dz  in kV/m — central difference in Excel
        if i == 1
            # forward difference: -(B3-B2)/((A3-A2)*1e-3) / 1e3
            #  z is in mm → Δz in m = Δz_mm * 1e-3
            sh2["C$r"] = XLSX.Formula("=-(B3-B2)/((A3-A2)*1E-3)/1000")
        elseif i == N_cells
            # backward difference
            prev = r - 1
            sh2["C$r"] = XLSX.Formula("=-(B$(r)-B$(r-1))/((A$(r)-A$(r-1))*1E-3)/1000")
        else
            nxt = r + 1
            prv = r - 1
            sh2["C$r"] = XLSX.Formula("=-(B$nxt-B$prv)/((A$nxt-A$prv)*1E-3)/1000")
        end
        sh2["D$r"] = E_field[i] / 1e3   # Julia cross-check value
    end

    # ════════════════════════════════════════════════════════════════════════
    # SHEET 3 — Ion_Velocity
    # Primitives  : z_mm (A), phi_V (B)
    # Formula cols:
    #   phi_max (C) = MAX(B:B)  — constant column, one formula per row
    #   dV_V    (D) = MAX(phi_max - phi_V, 0)
    #   v_ion_kms(E)= SQRT(2 * e_q * dV / Mi) / 1000
    # ════════════════════════════════════════════════════════════════════════
    sh3 = XLSX.addsheet!(xf, "Ion_Velocity")
    sh3["A1"] = "z_mm"
    sh3["B1"] = "phi_V"
    sh3["C1"] = "phi_max_V [formula]"
    sh3["D1"] = "dV_V [formula]"
    sh3["E1"] = "v_ion_kms [formula]"
    sh3["F1"] = "v_ion_kms [Julia ref]"
    sh3["G1"] = "Formula: v = SQRT(2*e_q*dV/Mi)/1000"

    for i in 1:N_cells
        r = i + 1
        sh3["A$r"] = z_mm[i]
        sh3["B$r"] = phi_V[i]
        # phi_max = MAX over entire B column (rows 2 to N+1)
        sh3["C$r"] = XLSX.Formula("=MAX(\$B\$2:\$B\$$(N_cells+1))")
        # dV = MAX(phi_max - phi(z), 0)
        sh3["D$r"] = XLSX.Formula("=MAX(C$r-B$r,0)")
        # v_ion [m/s] = sqrt(2 * 1.60218e-19 * dV / 2.18014e-25) / 1000
        sh3["E$r"] = XLSX.Formula("=SQRT(2*Constants!\$B\$2*D$r/Constants!\$B\$4)/1000")
        sh3["F$r"] = v_ion[i] / 1e3
    end

    # ════════════════════════════════════════════════════════════════════════
    # SHEET 4 — Ionisation_Source
    # Primitives  : z_mm (A), Te_eV (B), nn_m3 (C), ni_m3 (D)
    # Formula cols:
    #   k_iz (E) = Goebel & Katz (2008) App. E piecewise fit:
    #              IF(Te<5): 1e-20*(3.97+0.643*Te-0.0368*Te²)*exp(-12.127/Te)*v̄ₑ
    #              IF(Te≥5): 1e-20*(-1.031e-4*Te²+6.386*exp(-12.127/Te))*v̄ₑ
    #              v̄ₑ = SQRT(8*e_q*Te/(PI()*me))
    #   S_iz (F) = k_iz * nn * ni  [m⁻³s⁻¹]
    #   S_iz_1e22(G) = F / 1e22
    #   weight (H) = S_iz / SUM(S_iz)
    # ════════════════════════════════════════════════════════════════════════
    sh4 = XLSX.addsheet!(xf, "Ionisation_Source")
    sh4["A1"] = "z_mm"
    sh4["B1"] = "Te_eV"
    sh4["C1"] = "nn_m3"
    sh4["D1"] = "ni_m3"
    sh4["E1"] = "k_iz [formula]"
    sh4["F1"] = "S_iz [formula]"
    sh4["G1"] = "S_iz_1e22 [formula]"
    sh4["H1"] = "weight [formula]"
    sh4["I1"] = "k_iz [Julia ref]"
    sh4["J1"] = "S_iz [Julia ref]"
    sh4["K1"] = "Formula notes"
    sh4["K2"] = "k_iz = Goebel&Katz(2008) AppE piecewise IF(Te<5,...) formula  |  S_iz = k_iz*nn*ni  |  weight = S_iz/SUM(S_iz)"

    for i in 1:N_cells
        r = i + 1
        sh4["A$r"] = z_mm[i]
        sh4["B$r"] = Te_eV[i]
        sh4["C$r"] = nn_m3[i]
        sh4["D$r"] = ni_m3[i]
        # k_iz — Goebel & Katz (2008) Appendix E piecewise fit
        # Te floored at 0.5 eV (formula lower bound); v̄ₑ = SQRT(8*e_q*Te/(PI()*me))
        # Constants!$B$2 = e_q, Constants!$B$3 = me, Constants!$B$7 = Te threshold (5 eV)
        sh4["E$r"] = XLSX.Formula(
            "=IF(MAX(B$r,0.5)<Constants!\$B\$7," *
                "1E-20*(3.97+0.643*MAX(B$r,0.5)-0.0368*MAX(B$r,0.5)^2)*EXP(-Constants!\$B\$8/MAX(B$r,0.5))," *
                "1E-20*(-1.031E-4*MAX(B$r,0.5)^2+6.386*EXP(-Constants!\$B\$8/MAX(B$r,0.5))))" *
                "*SQRT(8*Constants!\$B\$2*MAX(B$r,0.5)/(PI()*Constants!\$B\$3))")
        # S_iz = k_iz * nn * ni
        sh4["F$r"] = XLSX.Formula("=E$r*C$r*D$r")
        # S_iz in units of 1e22
        sh4["G$r"] = XLSX.Formula("=F$r/1E22")
        # weight = S_iz / SUM(all S_iz)
        sh4["H$r"] = XLSX.Formula("=F$r/SUM(\$F\$2:\$F\$$(N_cells+1))")
        # Julia cross-check
        sh4["I$r"] = k_iz[i]
        sh4["J$r"] = S_iz[i]
    end

    # ════════════════════════════════════════════════════════════════════════
    # SHEET 5 — Hall_Parameter
    # Primitives  : z_B_mm (A), B_T (B), nu_en_interp (C)
    #   [nu_en_interp is simulation-grid-dependent; can't recompute in Excel
    #    without the full nn profile interpolation, so it stays a primitive]
    # Formula cols:
    #   omega_ce (D) = e_q * B / me
    #   Omega_e  (E) = omega_ce / (nu_en + 1)
    # ════════════════════════════════════════════════════════════════════════
    sh5 = XLSX.addsheet!(xf, "Hall_Parameter")
    sh5["A1"] = "z_B_mm"
    sh5["B1"] = "B_T"
    sh5["C1"] = "nu_en_interp_s"
    sh5["D1"] = "omega_ce [formula]"
    sh5["E1"] = "Omega_e [formula]"
    sh5["F1"] = "omega_ce [Julia ref]"
    sh5["G1"] = "Omega_e [Julia ref]"
    sh5["H1"] = "Formula notes"
    sh5["H2"] = "omega_ce = e_q*B/me  |  Omega_e = omega_ce/(nu_en+1)"

    for i in 1:N_B
        r = i + 1
        sh5["A$r"] = z_B_mm[i]
        sh5["B$r"] = B_T[i]
        sh5["C$r"] = nu_en_interp[i]
        # omega_ce = e * B / me
        sh5["D$r"] = XLSX.Formula("=Constants!\$B\$2*B$r/Constants!\$B\$3")
        # Omega_e = omega_ce / (nu_en + 1)
        sh5["E$r"] = XLSX.Formula("=D$r/(C$r+1)")
        sh5["F$r"] = omega_ce[i]
        sh5["G$r"] = Omega_e[i]
    end

    # ════════════════════════════════════════════════════════════════════════
    # SHEET 6 — Larmor_Radius
    # Primitives  : z_B_mm (A), B_T (B), Te_interp_eV (C)
    # Formula cols:
    #   Te_J     (D) = Te_eV * e_q
    #   r_L_m    (E) = SQRT(2*me*Te_J) / (e_q*B + 1e-10)
    #   r_L_mm   (F) = r_L_m * 1000
    #   r_L_ratio(G) = r_L_mm / channel_width_mm
    # ════════════════════════════════════════════════════════════════════════
    sh6 = XLSX.addsheet!(xf, "Larmor_Radius")
    sh6["A1"] = "z_B_mm"
    sh6["B1"] = "B_T"
    sh6["C1"] = "Te_interp_eV"
    sh6["D1"] = "Te_J [formula]"
    sh6["E1"] = "r_L_m [formula]"
    sh6["F1"] = "r_L_mm [formula]"
    sh6["G1"] = "r_L_ratio [formula]"
    sh6["H1"] = "r_L_mm [Julia ref]"
    sh6["I1"] = "Formula notes"
    sh6["I2"] = "r_L = SQRT(2*me*Te_J)/(e_q*B+1e-10)*1000  |  ratio = r_L_mm / ch_width_mm"

    for i in 1:N_B
        r = i + 1
        sh6["A$r"] = z_B_mm[i]
        sh6["B$r"] = B_T[i]
        sh6["C$r"] = Te_J_interp[i] / e_q   # store as eV (primitive)
        # Te in Joules
        sh6["D$r"] = XLSX.Formula("=C$r*Constants!\$B\$2")
        # r_L in metres
        sh6["E$r"] = XLSX.Formula("=SQRT(2*Constants!\$B\$3*D$r)/(Constants!\$B\$2*B$r+1E-10)")
        # r_L in mm
        sh6["F$r"] = XLSX.Formula("=E$r*1000")
        # r_L / channel_width
        sh6["G$r"] = XLSX.Formula("=F$r/Constants!\$B\$12")
        sh6["H$r"] = r_L_mm[i]
    end

    # ════════════════════════════════════════════════════════════════════════
    # SHEET 7 — ExB_Drift
    # Primitives  : z_mm (A), E_kVm (B), B_on_cells_G (C)
    # Formula cols:
    #   E_Vm     (D) = E_kVm * 1000
    #   B_T_cell (E) = B_on_cells_G / 10000
    #   v_ExB    (F) = ABS(E_Vm) / (B_T_cell + 1e-10) / 1000  [km/s]
    # ════════════════════════════════════════════════════════════════════════
    sh7 = XLSX.addsheet!(xf, "ExB_Drift")
    sh7["A1"] = "z_mm"
    sh7["B1"] = "E_kVm"
    sh7["C1"] = "B_on_cells_G"
    sh7["D1"] = "E_Vm [formula]"
    sh7["E1"] = "B_T_cell [formula]"
    sh7["F1"] = "v_ExB_kms [formula]"
    sh7["G1"] = "v_ExB_kms [Julia ref]"
    sh7["H1"] = "Formula notes"
    sh7["H2"] = "v_ExB = ABS(E[V/m]) / (B[T]+1e-10) / 1000  km/s"

    for i in 1:N_cells
        r = i + 1
        sh7["A$r"] = z_mm[i]
        sh7["B$r"] = E_field[i] / 1e3        # kV/m — primitive
        sh7["C$r"] = B_on_cells[i] * 1e4     # Gauss — primitive
        # Back to V/m
        sh7["D$r"] = XLSX.Formula("=B$r*1000")
        # Back to Tesla
        sh7["E$r"] = XLSX.Formula("=C$r/10000")
        # v_ExB [km/s]
        sh7["F$r"] = XLSX.Formula("=ABS(D$r)/(E$r+1E-10)/1000")
        sh7["G$r"] = v_ExB_kms[i]
    end

    # ════════════════════════════════════════════════════════════════════════
    # SHEET 8 — Scalar_Metrics  (all cross-sheet formulas — no hardcoded values)
    # Row layout: label | formula value | unit | formula string shown as text
    # ════════════════════════════════════════════════════════════════════════
    sh8 = XLSX.addsheet!(xf, "Scalar_Metrics")
    sh8["A1"] = "Parameter"
    sh8["B1"] = "Value [Excel formula]"
    sh8["C1"] = "Unit"
    sh8["D1"] = "Formula / Notes"

    # Helper: last data row in each sheet
    ecr = N_cells + 1   # last cell-grid row  (e.g. 53)
    ebr = N_B + 1       # last B-grid row     (e.g. 101)

    scalar_rows = [
        # ── Operating conditions ──────────────────────────────────────────
        ("── Operating conditions ──", "", "", ""),
        ("Discharge voltage Vd",
            "=Constants!B9",
            "V", "Constants!B9"),
        ("Mass flow rate mdot",
            "=Constants!B10*1E6",
            "mg/s", "Constants!B10 × 1e6"),
        ("Channel area A_ch",
            "=Constants!B11*1E4",
            "cm²", "Constants!B11 × 1e4"),
        ("", "", "", ""),

        # ── Potential & field ─────────────────────────────────────────────
        ("── Potential & field ──", "", "", ""),
        ("Anode potential φ_max",
            "=MAX(Electric_Field!B2:B$(ecr))",
            "V", "MAX(Electric_Field!phi_V)"),
        ("Exit potential φ_exit",
            "=Electric_Field!B$(ecr)",
            "V", "Electric_Field!phi_V[last row]"),
        ("Acceleration potential Δφ",
            "=B8-B9",
            "V", "φ_max − φ_exit  (B9−B10)"),
        ("Peak |E| E_max",
            "=MAX(ABS(Electric_Field!C2:C$(ecr)))*1000",
            "kV/m", "MAX(ABS(Electric_Field!E_kVm_formula))*1000"),
        ("z at E_max",
            "=INDEX(Electric_Field!A2:A$(ecr),MATCH(MAX(ABS(Electric_Field!C2:C$(ecr))),ABS(Electric_Field!C2:C$(ecr)),0))",
            "mm", "INDEX/MATCH on |E_kVm| column C"),
        ("", "", "", ""),

        # ── Ionisation zone ───────────────────────────────────────────────
        ("── Ionisation zone ──", "", "", ""),
        ("S_iz peak",
            "=MAX(Ionisation_Source!F2:F$(ecr))",
            "m⁻³s⁻¹", "MAX(Ionisation_Source!S_iz formula col F)"),
        ("z at S_iz peak",
            "=INDEX(Ionisation_Source!A2:A$(ecr),MATCH(MAX(Ionisation_Source!F2:F$(ecr)),Ionisation_Source!F2:F$(ecr),0))",
            "mm", "INDEX/MATCH on S_iz col F"),
        ("", "", "", ""),

        # ── Derived performance ───────────────────────────────────────────
        ("── Derived performance ──", "", "", ""),
        ("Ion exit velocity v_exit (S_iz-weighted)",
            "=SUMPRODUCT(Ionisation_Source!H2:H$(ecr),Ion_Velocity!E2:E$(ecr))",
            "km/s", "SUMPRODUCT(weight, v_ion_kms)"),
        ("Specific impulse Isp",
            "=B19*1000/Constants!B5",
            "s", "v_exit_ms / g0  → B19*1000/g0"),
        ("", "", "", ""),

        # ── Hall & magnetisation ──────────────────────────────────────────
        ("── Hall & magnetisation ──", "", "", ""),
        ("Peak Hall parameter Ω_e",
            "=MAX(Hall_Parameter!E2:E$(ebr))",
            "—", "MAX(Hall_Parameter!Omega_e formula col E)"),
        ("Min Larmor radius r_L",
            "=MIN(Larmor_Radius!F2:F$(ebr))",
            "mm", "MIN(Larmor_Radius!r_L_mm formula col F)"),
        ("Channel width",
            "=Constants!B12",
            "mm", "Constants!ch_w_mm"),
        ("", "", "", ""),

        # ── E×B drift ─────────────────────────────────────────────────────
        ("── E×B drift ──", "", "", ""),
        ("Peak v_ExB",
            "=MAX(ExB_Drift!F2:F$(ecr))",
            "km/s", "MAX(ExB_Drift!v_ExB_kms formula col F)"),
        ("z at peak v_ExB",
            "=INDEX(ExB_Drift!A2:A$(ecr),MATCH(MAX(ExB_Drift!F2:F$(ecr)),ExB_Drift!F2:F$(ecr),0))",
            "mm", "INDEX/MATCH on v_ExB col F"),
    ]

    for (i, row) in enumerate(scalar_rows)
        r = i + 1
        sh8["A$r"] = row[1]
        if row[2] != "" && startswith(row[2], "=")
            sh8["B$r"] = XLSX.Formula(row[2])
        else
            sh8["B$r"] = row[2]
        end
        sh8["C$r"] = row[3]
        sh8["D$r"] = row[4]
    end

end  # XLSX.openxlsx

println("Excel saved → $(abspath(xlsx_file))")
println("  Sheets: Constants | Electric_Field | Ion_Velocity | Ionisation_Source |")
println("          Hall_Parameter | Larmor_Radius | ExB_Drift | Scalar_Metrics")
println("  Formula columns labelled [formula]; primitive inputs labelled with no suffix.")
println("  Julia [ref] columns provided for cross-check — values should match.")