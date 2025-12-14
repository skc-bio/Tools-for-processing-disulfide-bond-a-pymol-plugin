# Tools-for-processing-disulfide-bond-a-pymol-plugin
I used Gemini 3 to write a plugin for handling the disulfide bond deficiency issue after Swiss-Model homology modeling.This plugin is written in Python and relies on the PyMOL software. 
# SS-Manager: PyMOL Disulfide Bond Repair Suite (v16.0)

**SS-Manager** is a comprehensive PyMOL plugin designed to solve topology and geometry issues related to disulfide bonds, specifically for **Homology Modeling** and **MD Simulation Preparation**.

It addresses three common pain points:
1.  **Missing Topology:** PDB files lacking `SSBOND` or `CONECT` records.
2.  **Geometric Distortion:** Modeled missing loops causing disulfide partners to be too far apart (> 10Å) for standard energy minimization.
3.  **Format Incompatibility:** Exported PDBs missing standard `SSBOND` headers required by GROMACS (`pdb2gmx`) or Amber (`tleap`).

---

## 📦 Features

* **Topology First:** Uses a CSV file as the "Source of Truth" to manage bond definitions.
* **Surgical Repair:** Includes `ss_snap`, a reproducible rigid-body translation tool to fix massive geometric gaps (> 5Å).
* **Diagnostics:** Instantly measures all bond distances and flags them as **OK** (Green), **Stretched** (Yellow), or **Broken** (Red).
* **MD Ready:** Automatically writes standard `SSBOND` records into the PDB header upon saving.
* **Topology Transfer:** Clones connectivity from a template structure to a target model.

---

## 📥 Installation

1.  Download `ss_manager.py` to your working directory.
2.  Open PyMOL.
3.  Run the script from the command line:
    ```pymol
    run ss_manager.py
    ```

---

## 🛠️ Command Reference

| Command | Description | Use Case |
| :--- | :--- | :--- |
| **`ss_export`** | Exports geometric bonds (< 3.2Å) to a CSV file. | Creating a baseline topology from a crystal structure. |
| **`ss_import`** | Reads a CSV and forces bond creation in PyMOL. | Visualizing the defined topology. |
| **`ss_check_dist`**| Reads a CSV and measures actual distances in the target object. | **Diagnosis:** Identifies "Broken" (Red) bonds. |
| **`ss_snap`** | Performs rigid-body translation to align a loop CYS to a core CYS. | **Repair:** Fixing gaps > 5Å (Requires Sculpting after). |
| **`ss_save_pdb`** | Saves the PDB and inserts `SSBOND` headers. | **Export:** Final step for GROMACS/Amber compatibility. |
| **`ss_transfer`** | Clones topology from Source to Target regardless of distance. | Quick cloning when sequence identity is 100%. |
| **`ss_compare`** | Compares Reference vs. Target bonds with Sequence Similarity. | Checking model fidelity against a template. |
| **`autobond_ss`** | Automatically bonds CYS pairs < 3.0Å. | Quick fix for simple structures. |

---

## 🧪 Standard Repair Protocol (The "Homology Fix")

Use this workflow when you have a perfect template (e.g., `apo`) and a broken model (e.g., `activate` with missing loops).

### Step 1: Define the Truth (CSV)
Extract the correct topology from the reference structure.
```pymol
load apo.pdb
ss_export apo, bonds.csv
Action: Open bonds.csv in a text editor. Add any missing bonds manually based on UniProt data. Save it.

### Step 2: Diagnose the Model
Load the broken model and check for "Red" flags.
Input Example:ss_check_dist protein, bonds.csv
Output Example: A:60 - B:70 | 18.08 Å | BROKEN (Red) -> Needs Repair!

### Step 3: Surgical Repair (Snap & Sculpt)
For every Broken (Red) bond, perform the following:
1.Snap (Align): Instantly move the loop to the core.
  # Usage: ss_snap [Anchor Atom], [Moving Atom], [Moving Scope/Loop]
  ss_snap /protein//A/60/SG, /protein//B/70/SG, /protein//B/60-80
2.Heal (Sculpt): Fix the backbone breakage caused by the snap.

### Step 4: Finalize and Export
Apply the topology visually and save the simulation-ready file.
Example:
# 1. Apply topology (Visual check)
ss_import protein, bonds.csv
# 2. Verify distances (Should be Green or Yellow now)
ss_check_dist protein, bonds.csv
# 3. Save with SSBOND headers
ss_save_pdb protein, bonds.csv, activate_final.pdb

⚠️ FAQ
Q: Why are some bonds Pink/Magenta? A: This indicates the distance is > 3.0Å.
3.0 - 4.5Å (Yellow status): Acceptable. Standard Energy Minimization (EM) in GROMACS/Amber will fix this automatically.
> 5.0Å (Red status): Dangerous. The forcefield might crash. Use ss_snap to fix these.
Q: Does this work with AmberTools? A: Yes. The ss_save_pdb command writes standard SSBOND records.

