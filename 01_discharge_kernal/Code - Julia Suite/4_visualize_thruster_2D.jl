# visualize_thruster_2D.jl
# ─────────────────────────────────────────────────────────────────────────────
# 2-D axisymmetric cross-section visualiser for HallThruster.jl output.
#
# Prerequisites
#   • My_HET.jl  has already been run (or included) in the same REPL session.
#   • 2D_axisymmetric_cross_section.jl has been run → `viz_data` is in scope.
#
# Usage (single REPL session)
#   julia> include("My_HET.jl")
#   julia> include("2D_axisymmetric_cross_section.jl")
#   julia> include("Visualise thruster 2d.jl")
#
# Output
#   MyHET_2D_viz.png  — saved in the current working directory
# ─────────────────────────────────────────────────────────────────────────────

using Plots
using Plots.PlotMeasures
gr()

# ── 0. Unpack viz_data ────────────────────────────────────────────────────────
r_i = viz_data.r_inner
r_o = viz_data.r_outer
L_ch = viz_data.L_ch
z_min = viz_data.z_min
z_max = viz_data.z_max

z_B = viz_data.z_B .* 1e3    # m → mm
B_prof = viz_data.B_prof .* 1e4    # T → Gauss

z_cells = vec(viz_data.z_cells) .* 1e3
Te = vec(viz_data.Te)
phi = vec(viz_data.phi)
ji = vec(viz_data.ji)
ni = vec(viz_data.ni) ./ 1e17
nn = vec(viz_data.nn) ./ 1e19

Vd = viz_data.Vd
mdot_mgs = viz_data.mdot_mgs

# Geometry in mm
r_i_mm = r_i * 1e3
r_o_mm = r_o * 1e3
L_ch_mm = L_ch * 1e3
z_max_mm = z_max * 1e3

# ── ZONE BOUNDARIES — derived from My_HET.jl geometry (auto-updates) ─────────
# These track L, z_peak, sigma defined in My_HET.jl
# Neutral Fill : 0  →  z_neutral_end  (60% of channel — before ionization ramps up)
# Ionisation   : z_neutral_end  →  z_peak  (B-field peak = ionization zone peak)
# Acceleration : z_peak  →  L_ch  (steepest potential drop, ions accelerate)
# Plume        : L_ch  →  domain end

z_neutral_end_mm = 0.60 * L_ch_mm          # 60% of channel length
z_ionise_end_mm = z_peak * 1e3            # B-field peak location (from My_HET.jl)
z_accel_end_mm = L_ch_mm                 # channel exit plane

println("Zone boundaries (auto-derived from geometry):")
@printf("  Neutral Fill : 0  →  %.1f mm\n", z_neutral_end_mm)
@printf("  Ionisation   : %.1f  →  %.1f mm\n", z_neutral_end_mm, z_ionise_end_mm)
@printf("  Acceleration : %.1f  →  %.1f mm (exit)\n", z_ionise_end_mm, z_accel_end_mm)
@printf("  Plume        : %.1f  →  %.1f mm\n", z_accel_end_mm, z_max_mm)

# ── PLUME DIVERGENCE ANGLE — voltage utilisation method (Option 3) ────────────
#
# Physical basis:
#   The beam divergence half-angle is set by the fraction of electron thermal
#   energy (Te_exit) at the exit plane that was NOT converted to directed axial
#   ion kinetic energy. Energy conservation gives:
#
#       α = arccos( sqrt(1 - Te_exit / Vd) )
#
#   This is geometry-independent and auto-scales with any operating point —
#   it uses only simulation outputs Te and Vd, so it works unchanged for any
#   channel size or discharge voltage you simulate.
#
#   IMPORTANT: This gives the IDEAL axial beam divergence only.
#   Experimental measurements will be larger (~2–3×) due to:
#     • Charge-exchange (CEX) ions populating large off-axis angles
#     • Near-field measurement distances (beam hasn't collimated)
#     • Charge-state mixing (multiply charged ions)
#
# All three methods printed to console for cross-checking; only Option 3 on plot.

exit_idx  = argmin(abs.(z_cells .- L_ch_mm))   # cell closest to exit plane
Te_exit   = Te[exit_idx]                         # eV — electron temperature at exit
phi_exit  = phi[exit_idx]                        # V  — residual plasma potential at exit

# Option 3 — voltage utilisation via Te_exit (PRIMARY — shown on plot)
α_volt_deg = Te_exit < Vd ?
    acosd(sqrt(max(0.0, 1.0 - Te_exit / Vd))) : 90.0

# Option 1 — geometric near-field (domain-dependent, reference only)
r_c_mm     = (r_i_mm + r_o_mm) / 2.0
z_ref_mm   = z_max_mm - L_ch_mm
α_geom_deg = atand(r_c_mm / z_ref_mm)

# Option 2 — ji decay to 50% in the plume (shown on plot as second cone)
#
# Physical basis:
#   Past the exit plane, ji rises to a plateau then slowly decays as the beam
#   expands. The axial distance at which ji first reaches 50% of its plateau
#   value (measured from the exit) gives a near-field beam half-width estimate.
#   Combined with the beam centroid radius r_c, the half-angle is:
#
#       α_ji = arctan( r_c / Δz_50 )
#
#   where Δz_50 = z(ji=50% plateau) - L_ch.
#   This is geometry-independent in the sense that it uses the actual current
#   density profile, not the domain size — but it does depend on how far the
#   plume has expanded within the domain.
#
# Implementation: search only in the plume region (z > L_ch_mm).
plume_mask   = z_cells .>= L_ch_mm
ji_plume     = ji[plume_mask]
z_plume      = z_cells[plume_mask]
ji_plateau   = maximum(ji_plume)                          # plateau value in plume
half_thresh  = 0.50 * ji_plateau                          # 50% threshold
# Find first index in plume where ji rises to 50% of plateau
half_idx_rel = findfirst(ji_plume .>= half_thresh)
if !isnothing(half_idx_rel) && half_idx_rel > 1
    z_ji_half_mm = z_plume[half_idx_rel]
    Δz_50        = z_ji_half_mm - L_ch_mm                # axial distance from exit
    α_ji_deg     = Δz_50 > 0.0 ? atand(r_c_mm / Δz_50) : 90.0
else
    # ji already at plateau at exit — beam fully accelerated, use full plume length
    z_ji_half_mm = L_ch_mm
    Δz_50        = z_ref_mm
    α_ji_deg     = atand(r_c_mm / z_ref_mm)
end

# Option 3 console label (phi-residual, kept for reference)
α_phi_deg  = phi_exit > 0.0 ? acosd(sqrt(max(0.0, 1.0 - phi_exit / Vd))) : 0.0

println("\nPlume divergence estimates:")
@printf("  Exit-plane Te          : %.2f eV  (%.1f%% of Vd)\n", Te_exit,  100*Te_exit/Vd)
@printf("  Exit-plane phi         : %.2f V   (%.1f%% of Vd)\n", phi_exit, 100*phi_exit/Vd)
@printf("  ji plateau (plume)     : %.2f A/m²\n", ji_plateau)
@printf("  z at 50%% ji plateau   : %.1f mm  (Δz=%.1f mm from exit)\n", z_ji_half_mm, Δz_50)
@printf("  ─────────────────────────────────────────────────────────────────\n")
@printf("  [PLOT] Volt-util  α    : %.1f °   arccos(sqrt(1 - Te_exit/Vd))  ← Option 3\n", α_volt_deg)
@printf("  [PLOT] ji-decay   α    : %.1f °   arctan(r_c / Δz_50)           ← Option 2\n", α_ji_deg)
@printf("  [REF]  Geometric  α    : %.1f °   arctan(r_c / z_plume)          ← Option 1 (domain-dependent)\n", α_geom_deg)
@printf("  [REF]  phi-resid  α    : %.1f °   arccos(sqrt(1 - phi_exit/Vd)) ← phi method\n", α_phi_deg)
@printf("  ─────────────────────────────────────────────────────────────────\n")
@printf("  NOTE: Experimental plume angles include CEX ions and near-field\n")
@printf("        effects — typically 2–3x larger than voltage-util value.\n")

# ── 1. Colours ────────────────────────────────────────────────────────────────
BG = RGB(0.08, 0.08, 0.12)
WALL = RGB(0.55, 0.55, 0.60)

C_Te = RGB(1.00, 0.55, 0.10)
C_phi = RGB(0.20, 0.85, 0.40)
C_ji = RGB(0.60, 0.30, 1.00)
C_ni = RGB(0.10, 0.75, 1.00)
C_nn = RGB(1.00, 0.30, 0.45)
C_B = RGB(0.50, 0.90, 1.00)

# ── 2. Helpers ────────────────────────────────────────────────────────────────
norm01(v) = (lo=minimum(v); hi=maximum(v); hi > lo ? (v .- lo) ./ (hi - lo) : zeros(length(v)))

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

function filled_plot!(p, xv, yv, col; lbl="", lw=2.0)
    alpha_fill = 0.15
    fc = RGBA(red(col), green(col), blue(col), alpha_fill)
    plot!(p, xv, zeros(length(xv)); color=fc, linewidth=0, label=false,
        fillrange=yv, fillcolor=fc, fillalpha=alpha_fill)
    plot!(p, xv, yv; color=col, linewidth=lw, label=lbl)
end

# ── add_zones! now uses derived boundaries ────────────────────────────────────
function add_zones!(p)
    # Ionisation zone: yellow
    vspan!(p, [z_neutral_end_mm, z_ionise_end_mm];
        alpha=0.07, color=:yellow, label=false)
    # Acceleration zone: deepskyblue
    vspan!(p, [z_ionise_end_mm, z_accel_end_mm];
        alpha=0.10, color=:deepskyblue, label=false)
end

# ── add_geo_axis! — secondary right Y-axis showing channel radial geometry ───
# Overlays a twinx() axis on panel `p` with:
#   • ylims fixed to the full radial span (-r_o_mm → +r_o_mm)
#   • four faint WALL-coloured hlines at ±r_o and ±r_i
#   • right-side tick labels at exactly those four positions (in mm)
#   • a dimmed "r (mm)" axis label
# This gives every physics panel a spatial geometry reference without
# touching or rescaling the primary (physics) Y-axis.
function add_geo_axis!(p)
    pt = twinx(p)
    plot!(pt;
        ylims             = (-r_o_mm, r_o_mm),
        yticks            = ([-r_o_mm, -r_i_mm, r_i_mm, r_o_mm],
                             ["-r_o=$(round(-r_o_mm,digits=1))",
                              "-r_i=$(round(-r_i_mm,digits=1))",
                              " r_i=$(round(r_i_mm,digits=1))",
                              " r_o=$(round(r_o_mm,digits=1))"]),
        yguide            = "r (mm)",
        guidefontcolor    = RGBA(red(WALL), green(WALL), blue(WALL), 0.7),
        guidefontsize     = 7,
        tickfontcolor     = RGBA(red(WALL), green(WALL), blue(WALL), 0.7),
        tickfontsize      = 6,
        foreground_color_axis = RGBA(red(WALL), green(WALL), blue(WALL), 0.4),
        grid              = false,
        minorticks        = false,
        background_color_inside = :transparent,
        label             = false,
    )
    # Four horizontal reference lines: outer wall (top/bot), inner wall (top/bot)
    hline!(pt, [ r_o_mm];
        color=WALL, linewidth=1.2, linestyle=:solid, alpha=0.35, label=false)
    hline!(pt, [-r_o_mm];
        color=WALL, linewidth=1.2, linestyle=:solid, alpha=0.35, label=false)
    hline!(pt, [ r_i_mm];
        color=WALL, linewidth=0.9, linestyle=:dash,  alpha=0.30, label=false)
    hline!(pt, [-r_i_mm];
        color=WALL, linewidth=0.9, linestyle=:dash,  alpha=0.30, label=false)
    return pt
end

# ── 3. Panel A — Geometry schematic ──────────────────────────────────────────
println("Building Panel A: geometry schematic...")

z_ch = collect(range(0.0, L_ch_mm, length=200))

# ── Centerline at 0 — classical cross-section convention ─────────────────────
r_i_top =  r_i_mm
r_o_top =  r_o_mm
r_i_bot = -r_i_mm
r_o_bot = -r_o_mm

y_pad      = r_i_mm * 0.3
z_Bpeak_mm = z_B[argmax(viz_data.B_prof .* 1e4)]

z_neutral_label_mm = z_neutral_end_mm / 2
z_ionise_label_mm  = (z_neutral_end_mm + z_ionise_end_mm) / 2
z_accel_label_mm   = (z_ionise_end_mm  + z_accel_end_mm)  / 2
z_plume_label_mm   = (z_accel_end_mm   + z_max_mm)        / 2

pA = plot(;
    xlims         = (z_min*1e3 - 1, z_max_mm + 1),
    ylims         = (-r_o_mm - y_pad, r_o_mm + y_pad),
    xlabel        = "Axial position  z (mm)",
    ylabel        = "Radius  r (mm)  [from centreline]",
    title         = "HET Cross-Section  |  Vd=$(Int(round(Vd))) V,  mdot=$(round(mdot_mgs,digits=4)) mg/s\nα_volt=$(round(α_volt_deg,digits=1))°  (volt-util - Orange line)      α_ji=$(round(α_ji_deg,digits=1))°  (ji-decay - Blue line)",
    titlefontsize = 10,
    legend        = :right,          # centred on right edge, inside the plume region
    legendfontsize = 7,
    _plt...)

# ── Centreline ────────────────────────────────────────────────────────────────
hline!(pA, [0.0]; color=:white, linestyle=:dot, linewidth=0.8, alpha=0.4, label=false)

# ── Channel walls — upper and lower ──────────────────────────────────────────
plot!(pA, z_ch, fill(r_o_top, length(z_ch)); color=WALL, linewidth=4, label="Wall (BN)")
plot!(pA, z_ch, fill(r_o_bot, length(z_ch)); color=WALL, linewidth=4, label=false)
plot!(pA, z_ch, fill(r_i_top, length(z_ch)); color=WALL, linewidth=4, label=false)
plot!(pA, z_ch, fill(r_i_bot, length(z_ch)); color=WALL, linewidth=4, label=false)

# ── Channel fill and inner core ───────────────────────────────────────────────
vspan!(pA, [0.0, L_ch_mm]; alpha=0.08, color=:cyan,  label=false)
hspan!(pA, [-r_i_mm, r_i_mm]; alpha=0.15, color=:gray, label=false)

# ── Zone shading ──────────────────────────────────────────────────────────────
add_zones!(pA)

# ── Reference lines ───────────────────────────────────────────────────────────
vspan!(pA, [L_ch_mm, z_max_mm]; alpha=0.05, color=:white, label=false)
vline!(pA, [L_ch_mm];    color=:white, linestyle=:dash,    linewidth=1.2, alpha=0.5, label="Exit plane")
vline!(pA, [z_Bpeak_mm]; color=C_B,   linestyle=:dashdot,  linewidth=1.0, alpha=0.8, label="B_max")

# ── Dimension annotations — stacked below inner wall, right side ──────────────
annot_x     = z_max_mm - 2.0        # right edge, inside domain
annot_y     = -r_i_mm - 1.5         # just below lower inner wall
line_gap    = r_i_mm * 0.18          # vertical spacing between lines

annotate!(pA, annot_x, annot_y,
          text("r_o = $(round(r_o_mm, digits=1)) mm", 7, :white, :right))
annotate!(pA, annot_x, annot_y - line_gap,
          text("w  = $(round(r_o_mm - r_i_mm, digits=1)) mm", 7, :yellow, :right))
annotate!(pA, annot_x, annot_y - 2*line_gap,
          text("r_i = $(round(r_i_mm, digits=1)) mm", 7, :white, :right))

# ── L_ch annotation — below lower outer wall ──────────────────────────────────
annotate!(pA, L_ch_mm/2, -r_o_mm - y_pad*0.6,
          text("L_ch = $(Int(round(L_ch_mm))) mm", 7, :white, :center))

# ── Zone labels — above upper outer wall ──────────────────────────────────────
annotate!(pA, z_neutral_label_mm, r_o_mm + y_pad*0.55, text("Neutral Fill", 7, :white,        :center))
annotate!(pA, z_ionise_label_mm,  r_o_mm + y_pad*0.55, text("Ionisation",   7, :yellow,       :center))
annotate!(pA, z_accel_label_mm,   r_o_mm + y_pad*0.55, text("Accel.",        6, :deepskyblue, :center))
annotate!(pA, z_plume_label_mm,   r_o_mm + y_pad*0.55, text("Plume",         7, :white,       :center))

# ── Plume divergence cones ────────────────────────────────────────────────────
# Both cones originate from beam centroid radius r_c at the exit plane.
# Option 3 (volt-util): ideal beam from energy conservation — tighter, physics-grounded
# Option 2 (ji-decay):  near-field observable from actual ji profile — broader

z_cone = [L_ch_mm, z_max_mm]

# Option 3 — voltage utilisation (orange dotted)
r3_top = [r_c_mm,  r_c_mm + z_ref_mm * tand(α_volt_deg)]
r3_bot = [-r_c_mm, -(r_c_mm + z_ref_mm * tand(α_volt_deg))]

plot!(pA, z_cone, r3_top;
    color=:orange, linewidth=1.2, linestyle=:dot, alpha=0.80, label="Volt-util cone")
plot!(pA, z_cone, r3_bot;
    color=:orange, linewidth=1.2, linestyle=:dot, alpha=0.80, label=false)

# Option 2 — ji 50% decay (cyan dotted) — typically broader than volt-util
r2_top = [r_c_mm,  r_c_mm + z_ref_mm * tand(α_ji_deg)]
r2_bot = [-r_c_mm, -(r_c_mm + z_ref_mm * tand(α_ji_deg))]

plot!(pA, z_cone, r2_top;
    color=:deepskyblue, linewidth=1.2, linestyle=:dot, alpha=0.70, label="ji-decay cone")
plot!(pA, z_cone, r2_bot;
    color=:deepskyblue, linewidth=1.2, linestyle=:dot, alpha=0.70, label=false)

# ── Cone angle labels — clamped inside ylims, placed at mid-plume x ───────────
y_lim_top   = r_o_mm + y_pad          # top of Panel A ylims
y_label_max = y_lim_top - 2.5         # max safe y for annotation text
z_mid_plume = L_ch_mm + 0.35 * z_ref_mm   # ~35% into plume — left of legend

# volt-util label: along its own cone line, clamped
r3_at_mid   = r_c_mm + 0.35 * z_ref_mm * tand(α_volt_deg)
y_v         = min(r3_at_mid + 2.0, y_label_max)
# annotate!(pA, z_mid_plume, y_v,
#     text("α=$(round(α_volt_deg,digits=1))° (volt-util)", 7, :orange, :left))

# ji-decay label: along its own cone line, clamped
r2_at_mid   = r_c_mm + 0.35 * z_ref_mm * tand(α_ji_deg)
y_j         = min(r2_at_mid + 2.0, y_label_max - 5.0)   # 5 mm below volt label if both clamped
# annotate!(pA, z_mid_plume, y_j,
#     text("α=$(round(α_ji_deg,digits=1))° (ji-decay)", 7, :deepskyblue, :left))


# ── 4. Panel B — Magnetic field ───────────────────────────────────────────────
println("Building Panel B: magnetic field...")

pB = plot(z_B, B_prof; color=C_B, linewidth=2.0,
    xlabel="z (mm)", ylabel="B (Gauss)",
    title="Magnetic Field  B(z)", titlefontsize=10,
    label="B(z)", _plt...)
add_zones!(pB)
vline!(pB, [L_ch_mm]; color=:white, linestyle=:dash, linewidth=1.0, alpha=0.5, label=false)
scatter!(pB, [z_Bpeak_mm], [maximum(B_prof)];
    color=C_B, markersize=5, markerstrokecolor=:white, label="B_peak")
add_geo_axis!(pB)

# ── 5. Panel C — Electron temperature ────────────────────────────────────────
println("Building Panel C: electron temperature...")

pC = plot(; xlabel="z (mm)", ylabel="Te (eV)",
    title="Electron Temperature  Te(z)", titlefontsize=10, _plt...)
filled_plot!(pC, z_cells, Te, C_Te; lbl="Te")
add_zones!(pC)
vline!(pC, [L_ch_mm]; color=:white, linestyle=:dash, linewidth=1.0, alpha=0.5, label=false)
z_Te_peak = z_cells[argmax(Te)]
scatter!(pC, [z_Te_peak], [maximum(Te)];
    color=C_Te, markersize=5, markerstrokecolor=:white, label="Peak")
annotate!(pC, z_Te_peak+1.5, maximum(Te)*0.93,
    text("$(round(maximum(Te),digits=1)) eV @ z=$(round(z_Te_peak,digits=1)) mm",
        7, :white, :left))
add_geo_axis!(pC)

# ── 6. Panel D — Plasma potential ────────────────────────────────────────────
println("Building Panel D: plasma potential...")

pD = plot(; xlabel="z (mm)", ylabel="phi (V)",
    title="Plasma Potential  phi(z)", titlefontsize=10, _plt...)
filled_plot!(pD, z_cells, phi, C_phi; lbl="phi")
add_zones!(pD)
vline!(pD, [L_ch_mm]; color=:white, linestyle=:dash, linewidth=1.0, alpha=0.5, label=false)
annotate!(pD, 2.0, maximum(phi)*0.92,
    text("$(round(maximum(phi),digits=0)) V -> $(round(minimum(phi),digits=0)) V",
        7, :white, :left))
add_geo_axis!(pD)

# ── 7. Panel E — Ion current density ─────────────────────────────────────────
println("Building Panel E: ion current density...")

pE = plot(z_cells, ji; color=C_ji, linewidth=2.0,
    xlabel="z (mm)", ylabel="ji (A/m2)",
    title="Ion Current Density  ji(z)", titlefontsize=10,
    label="ji", _plt...)
hline!(pE, [0.0]; color=:white, linestyle=:dot, linewidth=0.6, alpha=0.3, label=false)
add_zones!(pE)
vline!(pE, [L_ch_mm]; color=:white, linestyle=:dash, linewidth=1.0, alpha=0.5, label=false)
add_geo_axis!(pE)

# ── 8. Panel F — Ion number density ──────────────────────────────────────────
println("Building Panel F: ion number density...")

pF = plot(; xlabel="z (mm)", ylabel="ni (x1e17 m-3)",
    title="Ion Number Density  ni(z)", titlefontsize=10, _plt...)
filled_plot!(pF, z_cells, ni, C_ni; lbl="ni")
add_zones!(pF)
vline!(pF, [L_ch_mm]; color=:white, linestyle=:dash, linewidth=1.0, alpha=0.5, label=false)
add_geo_axis!(pF)

# ── 9. Panel G — Neutral number density ──────────────────────────────────────
println("Building Panel G: neutral number density...")

pG = plot(; xlabel="z (mm)", ylabel="nn (x1e19 m-3)",
    title="Neutral Density  nn(z)", titlefontsize=10, _plt...)
filled_plot!(pG, z_cells, nn, C_nn; lbl="nn")
add_zones!(pG)
vline!(pG, [L_ch_mm]; color=:white, linestyle=:dash, linewidth=1.0, alpha=0.5, label=false)
add_geo_axis!(pG)

# ── 10. Panel H — Normalised overlay ─────────────────────────────────────────
# No geo axis here — fixed 0–1 y-range makes radial mapping meaningless.
# The normalised overlay is already a cross-variable comparison panel.
println("Building Panel H: normalised overlay...")

pH = plot(; xlabel="z (mm)", ylabel="Normalised (a.u.)",
    title="All Profiles — Normalised Overlay", titlefontsize=10,
    ylims=(-0.05, 1.10), _plt...)
plot!(pH, z_cells, norm01(Te); color=C_Te, linewidth=1.8, label="Te")
plot!(pH, z_cells, norm01(phi); color=C_phi, linewidth=1.8, label="phi")
plot!(pH, z_cells, norm01(abs.(ji)); color=C_ji, linewidth=1.8, label="|ji|")
plot!(pH, z_cells, norm01(ni); color=C_ni, linewidth=1.8, label="ni")
plot!(pH, z_cells, norm01(nn); color=C_nn, linewidth=1.8, label="nn")
plot!(pH, z_B, norm01(B_prof); color=C_B, linewidth=1.5, linestyle=:dash, label="B")
add_zones!(pH)
vline!(pH, [L_ch_mm]; color=:white, linestyle=:dash, linewidth=1.0, alpha=0.5, label=false)

# ── 11. Compose and save ──────────────────────────────────────────────────────
println("Composing final figure...")

fig = plot(
    pA, pB, pC, pD, pE, pF, pG, pH;
    layout=(4, 2),
    size=(1600, 1800),
    background_color=BG,
    plot_title="MyHET_7mN  —  HallThruster.jl 2-D Axisymmetric Visualisation\n" *
        "Vd=$(Int(round(Vd))) V  |  mdot=$(round(mdot_mgs,digits=4)) mg/s  |  TwoZoneBohm(1/160,1/16)",
    plot_titlefontsize=11,
    plot_titlefontcolor=:white,
    left_margin=12mm,
    right_margin=14mm,   # bumped from 10 → 14 to give the geo-axis label room
    top_margin=10mm,
    bottom_margin=10mm,
    dpi=150,
)

outfile = "MyHET_2D_viz.png"
savefig(fig, outfile)
println("\nSaved: $(abspath(outfile))")
println("Zone summary:")
@printf("  Neutral Fill : 0 → %.1f mm  (%.0f%% of L_ch)\n", z_neutral_end_mm, 60.0)
@printf("  Ionisation   : %.1f → %.1f mm\n", z_neutral_end_mm, z_ionise_end_mm)
@printf("  Acceleration : %.1f → %.1f mm\n", z_ionise_end_mm, z_accel_end_mm)
display(fig)
