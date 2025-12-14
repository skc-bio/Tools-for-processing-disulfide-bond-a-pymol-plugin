# Tools-for-processing-disulfide-bond-a-pymol-plugin
I used Gemini 3 to write a plugin for handling the disulfide bond deficiency issue after Swiss-Model homology modeling.This plugin is written in Python and relies on the PyMOL software. 
# SS-Manager: PyMOL Disulfide-Bond Repair Suite (v16.0)
**SS-Manager** is a comprehensive PyMOL plugin that fixes topology and geometry problems around disulfide bonds—built for **homology modelling** and **MD-simulation prep**.
It solves three daily headaches:
1. **Missing topology** – PDB files that lack `SSBOND` or `CONECT` records  
2. **Geometric distortion** – modelled loops place the two CYS partners > 10 Å apart, beyond the reach of normal energy minimisation  
3. **Format incompatibility** – exported PDBs omit the standard `SSBOND` header that GROMACS (`pdb2gmx`) or Amber (`tleap`) insist on
---
## 📦 Features
* **Topology first** – a simple CSV file is the single source of truth  
* **Surgical repair** – `ss_snap` performs a reproducible rigid-body translation to close gaps > 5 Å  
* **Instant diagnostics** – bond lengths are colour-coded: **OK** (green), **Stretched** (yellow), **Broken** (red)  
* **MD-ready export** – `ss_save_pdb` automatically inserts correct `SSBOND` headers  
* **Topology transfer** – clone connectivity from a template to a target model in one step
---
## 📥 Installation
1. Download `ss_manager.py` to any folder  
2. Open PyMOL  
3. Run  
   ```pymol
   run ss_manager.py
   ```
---
## 🛠️ Command Reference
| Command | One-line description | Typical use |
|---|---|---|
| `ss_export` | Export bonds < 3.2 Å to CSV | Create a “truth” topology from a crystal structure |
| `ss_import` | Read CSV and force bonds in PyMOL | Visual check of the defined topology |
| `ss_check_dist` | Measure real distances vs CSV | Diagnose “red” bonds |
| `ss_snap` | Rigid-body translation to align loop CYS onto core CYS | Fix gaps > 5 Å (sculpt afterward) |
| `ss_save_pdb` | Save PDB with `SSBOND` headers | Final export for GROMACS / Amber |
| `ss_transfer` | Copy connectivity source→target regardless of distance | Quick clone when sequences are identical |
| `ss_compare` | Compare ref vs target bonds + sequence similarity | Validate model fidelity |
| `autobond_ss` | Auto-connect CYS pairs < 3.0 Å | One-click fix for simple cases |
---
## 🧪 Standard Repair Protocol (“Homology Fix”)
Use this when you have a perfect template (e.g. `apo.pdb`) and a broken model (e.g. `activate.pdb` with missing loops).
### Step 1 – Define the truth (CSV)
```pymol
load protein.pdb
ss_export protein, bonds.csv
```
Open `bonds.csv`, add any missing bonds from UniProt, save.
### Step 2 – Diagnose the model
```pymol
load activate.pdb
ss_check_dist activate, bonds.csv
```
Example output:  
`A:60 - B:70 | 18.08 Å | BROKEN (Red)` → repair needed
### Step 3 – Surgical repair (Snap & Sculpt)
For every **Red** bond:
1. **Snap** (align)
   ```pymol
   ss_snap /protein//A/60/SG, /protein//B/70/SG, /protein//B/60-80
   ```
2. **Sculpt** (heal) – run PyMOL sculpt or external EM to relax the backbone
### Step 4 – Finalise and export
```pymol
ss_import  protein, bonds.csv      # visual check
ss_check_dist protein, bonds.csv   # should be Green/Yellow
ss_save_pdb protein, bonds.csv, activate_final.pdb
```
---
## ⚠️ FAQ
**Q: Why are some bonds pink/magenta?**  
A: Distance > 3.0 Å  
- 3.0–4.5 Å (yellow): acceptable, standard EM will relax it  
- > 5 Å (red): dangerous, may crash the force-field—use `ss_snap`

**Q: Does this work with AmberTools?**  
A: Yes. `ss_save_pdb` writes standard `SSBOND` records that `tleap` recognises.
```
