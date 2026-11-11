# Code Documentation — 01 Discharge Kernel

This folder contains the markdown documentation for each Julia script in the
[Code - Julia Suite](../Code%20-%20Julia%20Suite/) folder.
Each `.md` file is a full walkthrough of its corresponding `.jl` script.

---

## Contents

| File | Covers | Description |
|---|---|---|
| `myHET_7mN.md` | `My_HET.jl` | Full code walkthrough — geometry, B-field, Config, SimParams, thrust sweep, 4-panel plots, Excel export, and design scaling reference |
| `2D_axisymmetric_cross_section.md` | `2D axisymmetric cross section.jl` | Data extraction walkthrough — reads plasma profiles from `sol_avg`, prints diagnostic summaries, and bundles all quantities into the `viz_data` named tuple for the downstream visualiser |

---

## How to Use

Each documentation file mirrors the section structure of its script.
Read the `.md` file alongside the `.jl` file to understand what each block does and why.

For the correct **run order** of the scripts, see the
[Code - Julia Suite README]([01_discharge_kernal/Code - Julia Suite](https://github.com/Shadow-Plasma/HET-Discharge-Solver/blob/93fd46e37863e9ca85d6524ceab90f6cf96a53e2/01_discharge_kernal/Code%20-%20Julia%20Suite/1_Code_Julia_Suite_README.md)).
