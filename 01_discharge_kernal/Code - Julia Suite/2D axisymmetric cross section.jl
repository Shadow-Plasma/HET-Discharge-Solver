# 2D axisymmetric cross section.jl
# No include! Variables are already in memory from Step 1 (My_HET.jl).
# Run My_HET.jl first in the same REPL session, then run this file.

using Printf

# ── 1. GEOMETRY ───────────────────────────────────────────────────────────────
r_inner = my_geometry.inner_radius    # m
r_outer = my_geometry.outer_radius    # m
L_ch = my_geometry.channel_length  # m
z_domain = my_config.domain            # (z_start, z_end) tuple

@printf("── Geometry ──────────────────────────────\n")
@printf("  Inner radius   : %.1f mm\n", r_inner * 1e3)
@printf("  Outer radius   : %.1f mm\n", r_outer * 1e3)
@printf("  Channel width  : %.1f mm\n", (r_outer - r_inner) * 1e3)
@printf("  Channel length : %.1f mm\n", L_ch * 1e3)
@printf("  Domain         : %.0f → %.0f mm\n", z_domain[1]*1e3, z_domain[2]*1e3)

# ── 2. MAGNETIC FIELD ─────────────────────────────────────────────────────────
z_B = z_vec                         # axial positions (m), length 100
B_prof = B_vec                         # B magnitude (T), length 100
B_peak = maximum(B_prof)
z_peak_B = z_B[argmax(B_prof)]

@printf("\n── Magnetic Field ────────────────────────\n")
@printf("  B_max          : %.1f Gauss  (%.4f T)\n", B_peak*1e4, B_peak)
@printf("  Peak location  : z = %.2f mm  (%.1f%% of channel)\n",
    z_peak_B*1e3, z_peak_B/L_ch*100)
@printf("  z_vec length   : %d points\n", length(z_B))

# ── 3. PLASMA PROFILE DATA (from sol_avg) ─────────────────────────────────────
z_cells = sol_avg.grid                 # axial cell centres (m)
Te_prof = sol_avg[:Tev][1]            # electron temperature (eV)
phi_prof = sol_avg[:potential][1]     # plasma potential (V)
ji_prof = sol_avg[:ji][1]           # ion current density (A/m²)

# Derived: neutral and ion number density if available
# (HallThruster.jl stores these as :nn and :ni)
ni_prof = sol_avg[:ni][1]           # ion number density (m⁻³)
nn_prof = sol_avg[:nn][1]           # neutral number density (m⁻³)

@printf("\n── Plasma Profile ────────────────────────\n")
@printf("  Grid cells     : %d points\n", length(z_cells))
@printf("  z range        : %.2f → %.2f mm\n", minimum(z_cells)*1e3, maximum(z_cells)*1e3)
@printf("  Te  range      : %.1f → %.1f eV  (peak at z = %.2f mm)\n",
    minimum(Te_prof), maximum(Te_prof),
    z_cells[argmax(Te_prof)]*1e3)
@printf("  φ (Plasma potential) range      : %.1f → %.1f V\n", minimum(phi_prof), maximum(phi_prof))
@printf("  ji (Ion Current) range      : %.2e → %.2e A/m²\n", minimum(ji_prof), maximum(ji_prof))
@printf("  ni (Ion density) range      : %.2e → %.2e m⁻³\n", minimum(ni_prof), maximum(ni_prof))
@printf("  nn (Neutral density) range      : %.2e → %.2e m⁻³\n", minimum(nn_prof), maximum(nn_prof))

# ── 4. OPERATING CONDITIONS ───────────────────────────────────────────────────
# In this HallThruster.jl version, mass flow rate is not a direct Config field.
# It lives in Config.propellants — try both; fall back to the known session value.
Vd = my_config.discharge_voltage         # V

mdot = try
    # v0.23+ stores propellant configs as a vector; first entry has the flow rate
    my_config.propellants[1].flow_rate
catch
    # Fallback: use the value from My_HET.jl (same REPL session, still in scope)
    0.6500071e-6
end

@printf("\n── Operating Conditions ──────────────────\n")
@printf("  Discharge voltage : %.0f V\n", Vd)
@printf("  Mass flow rate    : %.4f mg/s\n", mdot*1e6)

# ── 5. BUNDLE INTO A NAMED TUPLE FOR THE VISUALIZER ──────────────────────────
# Everything the 2D visualizer needs, in one clean struct
viz_data = (
    # Geometry
    r_inner=r_inner,
    r_outer=r_outer,
    L_ch=L_ch,
    z_min=z_domain[1],
    z_max=z_domain[2],

    # B-field
    z_B=z_B,
    B_prof=B_prof,

    # Plasma (axial)
    z_cells=z_cells,
    Te=Te_prof,
    phi=phi_prof,
    ji=ji_prof,
    ni=ni_prof,
    nn=nn_prof,

    # Operating conditions
    Vd=Vd,
    mdot_mgs=mdot * 1e6,
)

println("\n── viz_data bundle ready ─────────────────")
println("  Fields: ", keys(viz_data))
println("\nData extraction complete. Ready for visualizer.")