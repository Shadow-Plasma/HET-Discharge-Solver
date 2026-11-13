using HallThruster
using Printf
using Statistics

# ══════════════════════════════════════════════════════════════════════════════
# MyHET_7mN — CONFIG Generator V3
#
# PURPOSE
#   Given target T and Isp, find the (Vd, mdot) operating point and, if the
#   current geometry cannot hit the targets, adjust geometry within bounds.
#   Prints a validated CONFIG block to paste directly into My_HET.jl.
#
# SPEED DESIGN  (~30–60 s typical wall time)
#   • Bisection: 3 steps per geometry eval  (Isp vs Vd is smooth → fast bracket)
#   • Nelder-Mead: 25 outer evals max       (lean search, good enough for ±1%)
#   • Sim params: duration=2e-4, num_save=50 (half-lean → ~0.4–0.6 s per sim)
#   • Total sims: ≤ 25×3 + 7 = ~82         (well under 100)
#   • One final validation run at full lean fidelity to verify the result
#
# USAGE
#   1. Edit the TARGETS block below (T_target, Isp_target).
#   2. Paste your current My_HET.jl geometry into SEED GEOMETRY.
#   3. julia configgen_V3.jl
#   4. Copy the printed CONFIG block into My_HET.jl.
# ══════════════════════════════════════════════════════════════════════════════

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║                     TARGETS  —  EDIT HERE                              ║
# ╚══════════════════════════════════════════════════════════════════════════╝

const T_target = 7.30e-3   # target thrust           [N]
const Isp_target = 1065.0    # target specific impulse  [s]
const isp_tol = 8.0       # bisection tolerance      [s]  (wider → fewer sims)

# ── Objective weights ──────────────────────────────────────────────────────────
# obj = w_T·(ΔT/T)² + w_Isp·(ΔIsp/Isp)² − w_η·η_a + w_P·(ΔP/P_norm)²
const w_T = 2.0    # thrust match weight
const w_Isp = 2.0    # Isp match weight
const w_η = 0.5    # anode efficiency reward
const w_P = 0.3    # power penalty weight
const P_norm = 50.0   # power normalisation [W]

# ── Geometry search bounds ─────────────────────────────────────────────────────
#   Set lower==upper to FIX a parameter (removes it from optimisation).
const BOUNDS = Dict(
    :w_new => (0.007, 0.015),   # channel width         [m]   7–15 mm
    :ratio => (0.25, 0.45),    # w / R_mean
    :L_w_ratio => (1.8, 3.5),     # L_channel / w
    :B_ref => (0.012, 0.030),   # peak B at R_ref       [T]
    :z_peak_frac => (0.70, 0.95),    # B peak location / L
    :sigma_frac => (0.20, 0.45),    # B Gaussian σ / L
)

# ── Seed geometry — copy from your My_HET.jl CONFIG block ─────────────────────
const SEED_w_new = 0.010
const SEED_ratio = 0.31
const SEED_L_w_ratio = 2.5
const SEED_B_ref = 0.018
const SEED_z_peak_frac = 0.9
const SEED_sigma_frac = 0.3

# ── Fixed constants (match My_HET.jl exactly) ─────────────────────────────────
const R_ref = 0.0325    # reference radius for B scaling  [m]
const α_inner = 1/160
const α_outer = 1/16
const g₀ = 9.81      # ← matches My_HET.jl line 158 (uses 9.81, not 9.80665)

# ── Physical constants (for Larmor radius check) ──────────────────────────────
const e_charge = 1.602e-19   # elementary charge  [C]
const m_e = 9.109e-31   # electron mass      [kg]

# ── Vd search range ───────────────────────────────────────────────────────────
const Vd_min = 150.0
const Vd_max = 700.0

# ── Fast sim params (search phase) ────────────────────────────────────────────
#   duration=2e-4 + num_save=50 → ~0.4–0.6 s per sim on typical hardware
const FAST_dt = 1e-9
const FAST_duration = 2e-4
const FAST_num_save = 50

# ── Full lean sim params (validation run only) ────────────────────────────────
#   Matches My_HET.jl dt_lean / duration_lean / num_save_lean exactly
const FULL_dt = 1e-9
const FULL_duration = 4e-4
const FULL_num_save = 300

# ── Bisection steps (search phase) ────────────────────────────────────────────
#   3 steps = ±(Vd_range / 8) accuracy ≈ ±15–20 V → ±~15–20 s on Isp
#   Enough to steer NM; validation run is what actually counts.
const BISECT_STEPS_FAST = 3

# ── Nelder-Mead budget ────────────────────────────────────────────────────────
const NM_MAX_EVAL = 25     # 25 outer × 3 bisect = 75 fast sims ≈ 30–45 s
const NM_TOL = 5e-4   # stop early if simplex collapses

# ══════════════════════════════════════════════════════════════════════════════
#                       END OF USER SETTINGS
# ══════════════════════════════════════════════════════════════════════════════

# ── Parameter registry ────────────────────────────────────────────────────────
const PARAM_NAMES = [:w_new, :ratio, :L_w_ratio, :B_ref, :z_peak_frac, :sigma_frac]
const N_PARAMS = length(PARAM_NAMES)

# ── Unit-space helpers ────────────────────────────────────────────────────────
function to_phys(x::Vector{Float64})::Vector{Float64}
    [BOUNDS[n][1] + x[i] * (BOUNDS[n][2] - BOUNDS[n][1])
     for (i, n) in enumerate(PARAM_NAMES)]
end

function to_unit(phys::Vector{Float64})::Vector{Float64}
    [(phys[i] - BOUNDS[n][1]) / (BOUNDS[n][2] - BOUNDS[n][1])
     for (i, n) in enumerate(PARAM_NAMES)]
end

# ──────────────────────────────────────────────────────────────────────────────
# build_thruster
# ──────────────────────────────────────────────────────────────────────────────
function build_thruster(p::Vector{Float64})
    w_new, ratio, L_w_ratio, B_ref, z_peak_frac, sigma_frac = p

    R_mean = w_new / ratio
    geom = HallThruster.Geometry1D(
        inner_radius=R_mean - w_new/2,
        outer_radius=R_mean + w_new/2,
        channel_length=w_new * L_w_ratio,
    )
    L = geom.channel_length
    B_max = B_ref * (R_ref / R_mean)
    z_peak = z_peak_frac * L
    sigma = sigma_frac * L

    domain_end = 0.05
    z_vec = collect(range(0.0, domain_end, length=100))
    B_vec = B_max .* exp.(.-(z_vec .- z_peak) .^ 2 ./ (2*sigma^2))
    Bfield = HallThruster.MagneticField(file="", z=z_vec, B=B_vec)

    thr = HallThruster.Thruster(
        name="MyHET_7mN_Opt",
        geometry=geom,
        magnetic_field=Bfield,
        shielded=false,
    )
    return thr, geom
end

# ──────────────────────────────────────────────────────────────────────────────
# run_point: one simulation; fast=true for search phase, false for validation
# ──────────────────────────────────────────────────────────────────────────────
function run_point(thr, Vd::Float64, mdot::Float64; fast::Bool=true, return_sol::Bool=false)
    domain_end = 0.05
    cfg = HallThruster.Config(
        thruster=thr,
        discharge_voltage=Vd,
        anode_mass_flow_rate=mdot,
        domain=(0.0, domain_end),
        anom_model=HallThruster.TwoZoneBohm(α_inner, α_outer),
    )
    sim = HallThruster.SimParams(
        grid=HallThruster.EvenGrid(50),
        num_save=fast ? FAST_num_save : FULL_num_save,
        dt=fast ? FAST_dt : FULL_dt,
        duration=fast ? FAST_duration : FULL_duration,
    )
    sol = HallThruster.run_simulation(cfg, sim)
    if sol.retcode != :success
        return_sol && return (NaN, NaN, NaN, NaN, sol.retcode, nothing)
        return (NaN, NaN, NaN, NaN, sol.retcode)
    end
    avg = HallThruster.time_average(sol)
    T = HallThruster.thrust(avg)[]
    Id = HallThruster.discharge_current(avg)[]
    η = HallThruster.anode_eff(avg)[]
    Isp = T / (mdot * g₀)    # g₀=9.81, consistent with My_HET.jl
    P = Vd * Id
    return_sol && return (T, Isp, P, η, :success, sol)
    return (T, Isp, P, η, :success)
end

# ──────────────────────────────────────────────────────────────────────────────
# inner_bisect: bisect Vd to converge Isp to target.
# n_steps=BISECT_STEPS_FAST during search, more during validation.
# ──────────────────────────────────────────────────────────────────────────────
function inner_bisect(thr; fast::Bool=true, verbose::Bool=false)
    # Scaling-law seed: Isp ∝ √Vd, anchored at 235 V / 1061 s
    Vd_anchor = 235.0
    Isp_anchor = 1061.0
    Vd_est = clamp(Vd_anchor * (Isp_target / Isp_anchor)^2, Vd_min, Vd_max)
    mdot_est = T_target / (Isp_target * g₀)

    Vd_lo = max(Vd_min, Vd_est * 0.70)
    Vd_hi = min(Vd_max, Vd_est * 1.40)

    n_steps = fast ? BISECT_STEPS_FAST : 8
    Vd_best = Vd_est
    T_best = NaN
    Isp_best = NaN
    P_best = NaN
    η_best = NaN
    converged = false

    for iter in 1:n_steps
        Vd_mid = (Vd_lo + Vd_hi) / 2.0
        (T_sim, Isp_sim, P_sim, η_sim, rc) = run_point(thr, Vd_mid, mdot_est; fast=fast)

        if rc != :success || isnan(Isp_sim)
            verbose && @printf("    bisect %d: Vd=%.1f FAILED (%s)\n", iter, Vd_mid, rc)
            Vd_mid > Vd_best ? (Vd_hi = Vd_mid) : (Vd_lo = Vd_mid)
            continue
        end

        Δisp = Isp_sim - Isp_target
        if verbose
            @printf("    bisect %d: Vd=%5.1f V  T=%5.2f mN  Isp=%4.0f s  η=%.1f%%  ΔIsp=%+.0f s\n",
                iter, Vd_mid, T_sim*1e3, Isp_sim, η_sim*100, Δisp)
        end

        if isnan(Isp_best) || abs(Δisp) < abs(Isp_best - Isp_target)
            Vd_best = Vd_mid
            T_best = T_sim
            Isp_best = Isp_sim
            P_best = P_sim
            η_best = η_sim
        end

        if abs(Δisp) <= isp_tol; converged = true; break; end
        Isp_sim < Isp_target ? (Vd_lo = Vd_mid) : (Vd_hi = Vd_mid)
    end

    mdot_final = T_target / (Isp_best * g₀)
    return (Vd=Vd_best, mdot=mdot_final, T=T_best, Isp=Isp_best,
        P=P_best, η=η_best, converged=converged)
end

# ──────────────────────────────────────────────────────────────────────────────
# objective: score a unit-space point x ∈ [0,1]^N_PARAMS
# ──────────────────────────────────────────────────────────────────────────────
function objective(x::Vector{Float64}; verbose=false)
    x_c = clamp.(x, 0.0, 1.0)
    p = to_phys(x_c)
    thr, _ = try
        build_thruster(p)
    catch
        return (Inf, nothing)
    end

    r = inner_bisect(thr; fast=true, verbose=verbose)

    if any(isnan, [r.T, r.Isp, r.η, r.P])
        return (Inf, nothing)
    end

    ΔT = (r.T - T_target) / T_target
    ΔIsp = (r.Isp - Isp_target) / Isp_target
    ΔP = max(0.0, r.P - P_norm) / P_norm

    obj = w_T * ΔT^2 + w_Isp * ΔIsp^2 - w_η * r.η + w_P * ΔP^2
    return (obj, r)
end

# ══════════════════════════════════════════════════════════════════════════════
# Nelder-Mead in [0,1]^N_PARAMS
# BUG FIXED from V2: shrink loop now increments n_eval by 1 per vertex (not n).
# ══════════════════════════════════════════════════════════════════════════════
function nelder_mead(f_obj, x0::Vector{Float64};
    max_eval=NM_MAX_EVAL, tol=NM_TOL,
    α=1.0, γ=2.0, ρ=0.5, σ=0.5)
    n = length(x0)
    δ = 0.12   # 12% perturbation to build initial simplex

    simplex = [copy(x0) for _ in 1:(n+1)]
    for i in 1:n
        simplex[i+1][i] = clamp(x0[i] + δ, 0.0, 1.0)
    end

    vals = Vector{Float64}(undef, n+1)
    results = Vector{Any}(undef, n+1)
    for i in 1:(n+1)
        (vals[i], results[i]) = f_obj(simplex[i])
    end
    n_eval = n + 1

    best_obj = minimum(vals)
    no_improve_ctr = 0

    while n_eval < max_eval
        # 1. Order vertices best → worst
        order = sortperm(vals)
        simplex = simplex[order]
        vals = vals[order]
        results = results[order]

        # Early-stop if simplex has collapsed
        if abs(vals[end] - vals[1]) < tol
            no_improve_ctr += 1
            if no_improve_ctr >= 3; break; end
        else
            no_improve_ctr = 0
        end

        # 2. Centroid of all but worst
        x_o = mean(simplex[1:(end-1)])

        # 3. Reflection
        x_r = clamp.(x_o + α * (x_o - simplex[end]), 0.0, 1.0)
        (f_r, r_r) = f_obj(x_r)
        n_eval += 1

        if f_r < vals[1]
            # 4. Expansion
            x_e = clamp.(x_o + γ * (x_r - x_o), 0.0, 1.0)
            (f_e, r_e) = f_obj(x_e)
            n_eval += 1
            if f_e < f_r
                simplex[end] = x_e
                vals[end] = f_e
                results[end] = r_e
            else
                simplex[end] = x_r
                vals[end] = f_r
                results[end] = r_r
            end
        elseif f_r < vals[end-1]
            simplex[end] = x_r
            vals[end] = f_r
            results[end] = r_r
        else
            # 5. Contraction
            x_c = clamp.(x_o + ρ * (simplex[end] - x_o), 0.0, 1.0)
            (f_c, r_c) = f_obj(x_c)
            n_eval += 1
            if f_c < vals[end]
                simplex[end] = x_c
                vals[end] = f_c
                results[end] = r_c
            else
                # 6. Shrink  —  BUG FIX: n_eval += 1 per vertex (V2 had += n)
                for i in 2:(n+1)
                    simplex[i] = clamp.(simplex[1] + σ * (simplex[i] - simplex[1]), 0.0, 1.0)
                    (vals[i], results[i]) = f_obj(simplex[i])
                    n_eval += 1    # ← FIXED (was n_eval += n in V2)
                end
            end
        end

        new_best = minimum(vals)
        if new_best < best_obj - 1e-5
            best_obj = new_best
            bp = to_phys(clamp.(simplex[argmin(vals)], 0.0, 1.0))
            br = results[argmin(vals)]
            Vd_d = br !== nothing ? br.Vd : NaN
            η_d = br !== nothing ? br.η * 100 : NaN
            @printf("  [NM %3d]  obj=%.5f  Vd=%.1f V  η=%.1f%%  w=%.1fmm  B=%.3fT\n",
                n_eval, new_best, Vd_d, η_d, bp[1]*1000, bp[4])
        end
    end

    best_idx = argmin(vals)
    return (x_best=simplex[best_idx],
        obj_best=vals[best_idx],
        result=results[best_idx],
        n_eval=n_eval)
end

# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════
let
    t_start = time()

    println("\n" * repeat("═", 68))
    println("  MyHET_7mN  CONFIG Generator V3")
    println(repeat("═", 68))
    @printf("  Targets :  T = %.2f mN,  Isp = %.0f s\n", T_target*1e3, Isp_target)
    @printf("  g₀      :  %.2f m/s²  (matches My_HET.jl)\n", g₀)
    println()
    println("  Search space:")
    for n in PARAM_NAMES
        lo, hi = BOUNDS[n]
        if lo ≈ hi
            @printf("    %-14s  FIXED = %.4g\n", n, lo)
        else
            @printf("    %-14s  [%.4g, %.4g]\n", n, lo, hi)
        end
    end
    println()

    # ── PHASE 1: Evaluate seed ──────────────────────────────────────────────
    println(repeat("─", 68))
    println("  PHASE 1 — Seed evaluation (your current My_HET.jl geometry)")
    println(repeat("─", 68))

    x0 = to_unit([SEED_w_new, SEED_ratio, SEED_L_w_ratio,
        SEED_B_ref, SEED_z_peak_frac, SEED_sigma_frac])
    clamp!(x0, 0.0, 1.0)

    (seed_obj, seed_r) = objective(x0)
    if seed_r !== nothing
        @printf("  Seed:  obj=%.5f  T=%5.2f mN  Isp=%4.0f s  η=%4.1f%%  Vd=%5.1f V\n",
            seed_obj, seed_r.T*1e3, seed_r.Isp, seed_r.η*100, seed_r.Vd)

        ΔT_seed = abs(seed_r.T - T_target) / T_target * 100
        ΔIsp_seed = abs(seed_r.Isp - Isp_target) / Isp_target * 100
        if ΔT_seed < 2.0 && ΔIsp_seed < 1.0
            println("  → Seed geometry already meets targets. NM will refine operating point only.")
        else
            @printf("  → ΔT=%.1f%%  ΔIsp=%.1f%%  — geometry adjustment needed.\n",
                ΔT_seed, ΔIsp_seed)
        end
    else
        println("  Seed simulation failed — check bounds or geometry.")
    end

    # ── PHASE 2: Nelder-Mead ───────────────────────────────────────────────
    println()
    println(repeat("─", 68))
    @printf("  PHASE 2 — Nelder-Mead (max %d evals, fast sims)\n", NM_MAX_EVAL)
    println(repeat("─", 68))

    nm = nelder_mead(objective, x0)

    println()
    @printf("  NM done: %d evals,  best obj = %.5f\n", nm.n_eval, nm.obj_best)

    # ── PHASE 3: Full-fidelity validation ──────────────────────────────────
    println()
    println(repeat("─", 68))
    println("  PHASE 3 — Validation at full fidelity (duration=4e-4, num_save=300)")
    println(repeat("─", 68))

    p_best = to_phys(clamp.(nm.x_best, 0.0, 1.0))
    (thr_best, _geom) = build_thruster(p_best)
    r_val = inner_bisect(thr_best; fast=false, verbose=true)

    # ── Electron temperature & Larmor radius check ─────────────────────────
    # Re-run the best point at full fidelity to get raw sol for Te_peak
    w_val = p_best[1]   # channel width [m]
    B_ref_val = p_best[4]
    R_mean_val = w_val / p_best[2]
    B_peak_val = B_ref_val * (R_ref / R_mean_val)   # peak B at R_mean [T]

    Te_peak = NaN
    v_th = NaN
    r_larmor = NaN
    B_larmor_floor = NaN
    larmor_ok = false

    sol_val_tuple = run_point(thr_best, r_val.Vd, r_val.mdot; fast=false, return_sol=true)
    if length(sol_val_tuple) == 6 && sol_val_tuple[5] == :success
        sol_val = sol_val_tuple[6]
        avg_val = HallThruster.time_average(sol_val)
        Te_arr = avg_val[:Tev][1]          # electron temp profile [eV]
        Te_peak = maximum(Te_arr)           # peak Te [eV]
        v_th = sqrt(2 * e_charge * Te_peak / m_e)   # thermal speed [m/s]
        r_larmor = m_e * v_th / (e_charge * B_peak_val) # Larmor radius [m]
        # Design rule: r_larmor < w/10  (electrons confined, not streaming across)
        B_larmor_floor = 10 * m_e * v_th / (e_charge * w_val)  # min B for confinement [T]
        larmor_ok = r_larmor < w_val / 10
    end

    t_total = time() - t_start

    # ── Print optimised CONFIG block ───────────────────────────────────────
    w_new_opt, ratio_opt, L_w_ratio_opt,
    B_ref_opt, z_peak_frac_opt, sigma_frac_opt = p_best

    R_mean_opt = w_new_opt / ratio_opt
    L_opt = w_new_opt * L_w_ratio_opt

    println()
    println(repeat("═", 68))
    println("  OPTIMISED CONFIG — paste into My_HET.jl  (CONFIGURATION block)")
    println(repeat("═", 68))
    println()
    println("  # ── Geometry ─────────────────────────────────────────────────────")
    @printf("  const w_new     = %.4f   # channel width [m]  (%.1f mm)\n",
        w_new_opt, w_new_opt*1000)
    @printf("  const ratio     = %.4f   # w / R_mean\n", ratio_opt)
    @printf("  const L_w_ratio = %.3f   # L_channel / w\n", L_w_ratio_opt)
    println()
    println("  # ── Magnetic field ───────────────────────────────────────────────")
    @printf("  const B_ref       = %.5f  # peak B at R_ref [T]\n", B_ref_opt)
    @printf("  const R_ref       = %.4f  # fixed — do not change\n", R_ref)
    @printf("  const z_peak_frac = %.3f  # B peak / L\n", z_peak_frac_opt)
    @printf("  const sigma_frac  = %.3f  # B Gaussian σ / L\n", sigma_frac_opt)
    println()
    println("  # ── Nominal operating point ──────────────────────────────────────")
    @printf("  const Vd_nom   = %.1f     # discharge voltage [V]\n", r_val.Vd)
    @printf("  const mdot_ref = %.7e # anode mass flow [kg/s]  (%.4f mg/s)\n",
        r_val.mdot, r_val.mdot*1e6)
    println()
    println("  # ── Simulation-verified results ──────────────────────────────────")
    @printf("  # Thrust   : %6.2f mN   (target: %.2f mN,  error: %+.2f mN)\n",
        r_val.T*1e3, T_target*1e3, (r_val.T - T_target)*1e3)
    @printf("  # Isp      : %6.0f s    (target: %.0f s,   error: %+.0f s)\n",
        r_val.Isp, Isp_target, r_val.Isp - Isp_target)
    @printf("  # η_anode  : %6.1f %%\n", r_val.η*100)
    @printf("  # P_disch  : %6.1f W    (discharge power)\n", r_val.P)
    @printf("  # P_sys    : %6.1f – %.1f W   (η_PPU = 85–90%%)\n",
        r_val.P/0.90, r_val.P/0.85)
    println()
    println("  # ── Derived geometry (for reference) ────────────────────────────")
    @printf("  # R_mean   : %.2f mm\n", R_mean_opt*1000)
    @printf("  # R_inner  : %.2f mm\n", (R_mean_opt - w_new_opt/2)*1000)
    @printf("  # R_outer  : %.2f mm\n", (R_mean_opt + w_new_opt/2)*1000)
    @printf("  # L_channel: %.2f mm\n", L_opt*1000)
    @printf("  # w/R_mean : %.3f  (healthy: 0.25–0.35)\n", w_new_opt/R_mean_opt)
    @printf("  # L/w      : %.2f   (healthy: 2.0–3.0)\n", L_opt/w_new_opt)
    @printf("  # B_peak   : %.1f Gauss  (at R_mean)\n",
        B_ref_opt * (R_ref / R_mean_opt) * 1e4)
    println()
    println("  # ── Geometry change from seed ───────────────────────────────────")
    for (name, seed, opt, unit, scale) in [
        ("w_new", SEED_w_new, w_new_opt, "mm", 1000.0),
        ("ratio", SEED_ratio, ratio_opt, "—", 1.0),
        ("L_w_ratio", SEED_L_w_ratio, L_w_ratio_opt, "—", 1.0),
        ("B_ref", SEED_B_ref, B_ref_opt, "T", 1.0),
        ("z_peak_frac", SEED_z_peak_frac, z_peak_frac_opt, "—", 1.0),
        ("sigma_frac", SEED_sigma_frac, sigma_frac_opt, "—", 1.0),
    ]
        Δ = (opt - seed) * scale
        pct = abs(seed) > 1e-12 ? (opt - seed)/seed*100 : 0.0
        marker = abs(pct) < 0.5 ? "  (unchanged)" : @sprintf("  Δ=%+.3f %s (%+.1f%%)", Δ, unit, pct)
        @printf("  #   %-14s  seed=%.4g → opt=%.4g%s\n", name, seed*scale, opt*scale, marker)
    end
    println()
    println("  # ── Convergence ─────────────────────────────────────────────────")
    if r_val.converged
        println("  # [OK] Validation bisection converged within ±$(isp_tol) s")
    else
        @printf("  # [!]  Validation bisection did not converge: Isp=%.0f s (Δ=%+.0f s)\n",
            r_val.Isp, r_val.Isp - Isp_target)
        println("  #      If Δ < ±20 s the result is still usable.")
        println("  #      To improve: widen Vd_min/Vd_max or increase isp_tol.")
    end
    println()
    println("  # ── Optimiser stats ─────────────────────────────────────────────")
    @printf("  # NM evals : %d  (fast sims)\n", nm.n_eval)
    @printf("  # Obj value: %.5f  (0 = perfect T+Isp match + max η)\n", nm.obj_best)
    @printf("  # Wall time: %.0f s\n", t_total)
    println()
    println("  # ── Electron confinement check (Larmor radius) ──────────────────")
    if !isnan(Te_peak)
        @printf("  # Te_peak      : %.1f eV\n", Te_peak)
        @printf("  # v_th (elec)  : %.3e m/s\n", v_th)
        @printf("  # B_peak       : %.1f Gauss\n", B_peak_val * 1e4)
        @printf("  # r_Larmor     : %.3f mm\n", r_larmor * 1e3)
        @printf("  # w / 10       : %.3f mm  (confinement limit)\n", w_val * 1e3 / 10)
        @printf("  # B_min needed : %.1f Gauss  (for r_L < w/10)\n", B_larmor_floor * 1e4)
        if larmor_ok
            println("  # [OK] Electrons well-magnetised  (r_L < w/10)")
        else
            println("  # [!]  r_L ≥ w/10 — electrons under-magnetised; consider increasing B_ref")
        end
    else
        println("  # [--] Larmor check skipped (sol unavailable)")
    end
    println(repeat("═", 68))

end
