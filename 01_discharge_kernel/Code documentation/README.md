# Code Documentation — 01 Discharge Kernel

This folder contains the markdown documentation for each Julia script in the
[Code - Julia Suite](../Code%20-%20Julia%20Suite/) folder.
Each `.md` file is a full walkthrough of its corresponding `.jl` script.

---

## Contents

| # | File | Covers | Description |
|---|---|---|---|
| 0 | `0_ConfigV3.md` | `0_ConfigV3.jl` | Automated CONFIG generator — Nelder-Mead + bisection search for optimal (Vd, mdot, geometry). Run standalone before the main pipeline. |
| 1 | `1_myHET_7mN.md` | `1_My_HET.jl` | Full code walkthrough — geometry, B-field, Config, SimParams, thrust sweep, 4-panel plots, Excel export, and design scaling reference |
| 2 | `2_2D_axisymmetric_cross_section.md` | `2_2D axisymmetric cross section.jl` | Data extraction walkthrough — reads plasma profiles from `sol_avg`, prints diagnostic summaries, and bundles all quantities into the `viz_data` named tuple for the downstream visualiser |
| 3 | `3_Visualise_thruster_2d.md` | `3_Visualise thruster 2d.jl` | 8-panel 2D visualiser walkthrough — zone shading, plume divergence angles, geometry schematic, panel-by-panel breakdown |
| 4 | `4_Derived_Quantities.md` | `4_Derived_Quantities.jl` | Derived plasma physics walkthrough — E-field, ion velocity, Hall parameter, E×B drift, ionisation source, Excel export |

---

## How to Use

Each documentation file mirrors the section structure of its corresponding script.
Read the `.md` file alongside the `.jl` file to understand what each block does and why.

For the correct **run order** of the scripts, see the
[Code - Julia Suite README](../Code%20-%20Julia%20Suite/1_Code_Julia_Suite_README.md).
