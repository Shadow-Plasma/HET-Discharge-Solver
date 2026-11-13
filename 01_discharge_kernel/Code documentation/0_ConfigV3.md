# ConfigV3 — Automated CONFIG Generator

Given a thrust and Isp target, finds the optimal `(Vd, mdot)` operating point and, if the current geometry cannot hit the targets, adjusts geometry within user-defined bounds. Prints a validated CONFIG block ready to paste directly into `My_HET.jl`.

> **Prerequisites:** Run as a standalone script — does **not** require `My_HET.jl` to be in scope first. All dependencies are self-contained.

---

## Usage

```julia
# Edit the TARGETS block and SEED GEOMETRY inside the script, then run:
julia> include("ConfigV3.jl")
```

Copy the printed CONFIG block into the `# ── CONFIGURATION block` section of `My_HET.jl`.

---

## User Settings

All user-editable parameters are at the top of the script.

### Targets

```julia
const T_target   = 7.30e-3   # target thrust            [N]
const Isp_target = 1065.0    # target specific impulse  [s]
const isp_tol    = 8.0       # bisection Isp tolerance  [s]
```

### Objective Function

```
obj = w_T·(ΔT/T)² + w_Isp·(ΔIsp/Isp)² − w_η·η_a + w_P·(ΔP/P_norm)²
```

| Weight | Default | Purpose |
|---|---|---|
| `w_T` | 2.0 | Thrust match |
| `w_Isp` | 2.0 | Isp match |
| `w_η` | 0.5 | Anode efficiency reward |
| `w_P` | 0.3 | Power penalty |
| `P_norm` | 50 W | Power normalisation reference |

### Geometry Search Bounds

| Parameter | Range | Description |
|---|---|---|
| `w_new` | 7–15 mm | Channel width |
| `ratio` | 0.25–0.45 | `w / R_mean` |
| `L_w_ratio` | 1.8–3.5 | `L_channel / w` |
| `B_ref` | 0.012–0.030 T | Peak B at reference radius |
| `z_peak_frac` | 0.70–0.95 | B peak location as fraction of L |
| `sigma_frac` | 0.20–0.45 | B Gaussian σ as fraction of L |

> Set `lower == upper` to **fix** a parameter and remove it from the optimisation.

### Seed Geometry

Paste your current `My_HET.jl` geometry values here before running:

```julia
const SEED_w_new       = 0.010
const SEED_ratio       = 0.31
const SEED_L_w_ratio   = 2.5
const SEED_B_ref       = 0.018
const SEED_z_peak_frac = 0.9
const SEED_sigma_frac  = 0.3
```

---

## Algorithm

Three-phase search:

```
Phase 1 — Seed evaluation
          Evaluate the current My_HET.jl geometry at the target
          operating point. Reports ΔT% and ΔIsp% to decide if
          geometry adjustment is needed.
          ↓
Phase 2 — Nelder-Mead geometry search (fast sims)
          Optimises all 6 geometry parameters simultaneously.
          Each NM eval calls inner bisection to find the best
          Vd for that geometry.
          ↓
Phase 3 — Full-fidelity validation
          Re-runs the best geometry at full sim params to verify
          the result. Includes electron confinement (Larmor radius)
          check.
```

### Inner Bisection (Vd solver)

For a given geometry, bisects discharge voltage `Vd ∈ [Vd_min, Vd_max]` to converge `Isp → Isp_target`. Mass flow rate is set analytically:

```
mdot = T_target / (Isp_target × g₀)
```

Vd seed from scaling law:

```
Vd_est = Vd_anchor × (Isp_target / Isp_anchor)²

Vd_anchor = 235 V,  Isp_anchor = 1061 s
```

### Nelder-Mead

Operates in unit space `[0, 1]^6` (each parameter normalised to its bounds). Parameters:

| Parameter | Value |
|---|---|
| Max evals | 25 (fast sims) |
| Simplex perturbation δ | 12% |
| Stop tolerance | 5 × 10⁻⁴ |
| Reflection α | 1.0 |
| Expansion γ | 2.0 |
| Contraction ρ | 0.5 |
| Shrink σ | 0.5 |

> **V3 bug fix:** The shrink step now increments `n_eval` by 1 per vertex (V2 incorrectly incremented by `n`).

---

## Simulation Parameters

| Phase | `dt` | `duration` | `num_save` | Approx. time/sim |
|---|---|---|---|---|
| Fast (search) | 1 ns | 0.2 ms | 50 | ~0.4–0.6 s |
| Full (validation) | 1 ns | 0.4 ms | 300 | ~1–2 s |

Total simulation budget: ≤ 25 × 3 + 7 ≈ 82 fast sims → **~30–60 s typical wall time**.

---

## Electron Confinement Check

After the validation run, the script checks whether electrons are well-magnetised using the Larmor radius criterion:

```
r_L = m_e × v_th / (e × B_peak)

v_th = sqrt( 2 × e × Te_peak / m_e )

Design rule:  r_L < w / 10
```

If `r_L ≥ w/10` the script prints a warning and reports the minimum B needed for confinement.

---

## Console Output

```
════════════════════════════════════════════════════════════════════
  MyHET_7mN  CONFIG Generator V3
════════════════════════════════════════════════════════════════════
  Targets :  T = 7.30 mN,  Isp = 1065 s
  g₀      :  9.81 m/s²  (matches My_HET.jl)

  PHASE 1 — Seed evaluation
  ─────────────────────────────────────────────────────────────────
  Seed:  obj=X.XXXXX  T= X.XX mN  Isp=XXXX s  η=XX.X%  Vd=XXX.X V

  PHASE 2 — Nelder-Mead (max 25 evals, fast sims)
  ─────────────────────────────────────────────────────────────────
  [NM  XX]  obj=X.XXXXX  Vd=XXX.X V  η=XX.X%  w=XX.Xmm  B=X.XXX T

  PHASE 3 — Validation at full fidelity
  ─────────────────────────────────────────────────────────────────
  [bisection verbose output]

════════════════════════════════════════════════════════════════════
  OPTIMISED CONFIG — paste into My_HET.jl
════════════════════════════════════════════════════════════════════

  # Thrust   :   X.XX mN   (target: X.XX mN,  error: +X.XX mN)
  # Isp      :   XXXX s    (target: XXXX s,   error: +X s)
  # η_anode  :   XX.X %
  # P_disch  :   XX.X W
  # Wall time:   XX s
════════════════════════════════════════════════════════════════════
```

---

## Output

No files are written. The optimised CONFIG block is printed to the console for copy-paste into `My_HET.jl`.

---

## Physical Constants

| Constant | Symbol | Value | Unit |
|---|---|---|---|
| Elementary charge | `e_charge` | 1.602 × 10⁻¹⁹ | C |
| Electron mass | `m_e` | 9.109 × 10⁻³¹ | kg |
| Standard gravity | `g₀` | 9.81 | m/s² |
| Reference radius | `R_ref` | 32.5 mm | m |

> `g₀ = 9.81` matches `My_HET.jl` exactly. The standard value is 9.80665 m/s² — using 9.80665 here would cause a small but consistent Isp offset between the two scripts.

---

## Dependencies

```julia
using HallThruster   # v0.23.2
using Printf
using Statistics
```

No additional packages required beyond those used in `My_HET.jl`.

---

## Version History

| Version | Change |
|---|---|
| V3 | Bug fix: shrink step increments `n_eval` by 1 per vertex (V2 had `+= n`) |
| V2 | Added Larmor radius check; geometry search bounds |
| V1 | Initial bisection-only operating point finder |
