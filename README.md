# Lattice Dispersion Visualization (MATLAB)
This repository provides a MATLAB script for importing, analyzing, and visualizing the dynamic behavior of 2D rectangular lattice structures from Abaqus eigenfrequency analysis.

## Generating the `.rpt` files (Abaqus or COMSOL)

The scripts expect **plain-text, tab/space-delimited** reports with the following columns per k-point and mode (order can vary as long as columns are present):  
`kx, ky, mode, frequency, eff_mass_Tx, eff_mass_Ty, eff_mass_Rz`

**Units (recommended):** frequency in **Hz** (or rad/s if you keep it consistent everywhere) and wave vectors **kx, ky** in **1/m**. Use a BZ **PATH** (Γ–X–M–Γ) for `abaqus_Frequencies_PATH.rpt` and a uniform **GRID** of (kx, ky) for `abaqus_Frequencies_GRID.rpt`.

**Abaqus/CAE (Frequency step)**
1. Create a **Frequency** step for your unit cell. (If using periodic media, apply Bloch/Floquet-type periodic constraints and solve one job per k-point along the path or grid.)
2. Run the job(s). In **Visualization**:
   - `Report → Modal and participation factors…`
   - Select **All modes** and enable **Frequencies** and **Effective mass**; choose components **TX, TY** and **RZ**.
   - Set output to **tab-delimited text** and **Append** when collecting multiple k-points.
   - Save as `example_data/abaqus_Frequencies_PATH.rpt` (for the path sweep) or `example_data/abaqus_Frequencies_GRID.rpt` (for the grid sweep).
3. If your report tool doesn’t add `kx, ky` automatically, keep one file per k-point and concatenate them in order (PATH) or row-major (GRID), or include a short header line per block indicating `kx, ky`.

**COMSOL Multiphysics (Eigenfrequency + Bloch–Floquet)**
1. Solid Mechanics → add **Periodic Condition** with **Floquet (Bloch–Floquet)** phase shift; define the wave vector parameters **kx, ky**.
2. Study: **Eigenfrequency**; add a **Parametric Sweep** over the desired k-points (PATH or GRID). Choose number of modes and (if needed) a frequency shift.
3. After solving:
   - `Results → Derived Values → Global Evaluation`: add **eigfreq** (eigenfrequencies) and evaluate for all parameter sets.
   - `Results → Export → Data`: export the evaluation table as **tab-delimited text** and save as `example_data/abaqus_Frequencies_PATH.rpt` or `...GRID.rpt`.
   - *(Optional, for polarization plots)* Add global evaluations for mass-weighted displacement components to produce **eff_mass_Tx**, **eff_mass_Ty**, and an out-of-plane rotational measure (**eff_mass_Rz**) compatible with your Abaqus columns. If you prefer, you can export eigenvectors and compute these quantities offline to match the expected columns.

---

## Main features:
- **Path Analysis:** Frequency dispersion curves along high-symmetry paths in the Brillouin Zone (using `abaqus_Frequencies_PATH.rpt`).
- **Grid Analysis:** 3D dispersion surfaces over a grid of wave vectors (kx, ky) (using `abaqus_Frequencies_GRID.rpt`).
- **Polarization Visualization:** Effective mass data for translation & rotation modes.
- **Group Velocity Analysis:** Direction and speed of energy propagation.

## Getting Started
**Prepare Abaqus Output Files:**  
Export eigenfrequency data (including effective mass and rotational mass) as `.rpt` files.  
Required files:
- `abaqus_Frequencies_PATH.rpt` (for boundary/path analysis)
- `abaqus_Frequencies_GRID.rpt` (for full IBZ grid analysis)  
Place these files in the `example_data/` folder.

**Edit File Paths (if needed):**  
In the script, set the paths to your local `.rpt` files if different.

**Run the Script:**  
Open `Rectangular_lattice_3by1.m` in MATLAB and run it.

**Customize & Extend:**  
All code is commented for easy extension (e.g., different lattice geometries, advanced plots).

## Example Data
Example Abaqus output files are provided in `example_data/`:
- `abaqus_Frequencies_PATH.rpt` — for path analysis
- `abaqus_Frequencies_GRID.rpt` — for grid analysis

Use these files to test and reproduce the MATLAB analysis.

**Note:**
The code and `.rpt` files provided correspond to the rectangular lattice illustrated in **Figure 6** of the journal article "Identifying elastic wave polarization and bandgaps in periodic solid media" and to **Section 5.3.1** of the dissertation "A Novel Method to Study the Dynamic Behavior of Lattice Structures".

## Further reading
For additional background and methodology details, see the dissertation:

Carrillo-Munoz, Maria J. "A Novel Method to Study the Dynamic Behavior of Lattice Structures." PhD Dissertation, Wichita State University, 2023. https://soar.wichita.edu/server/api/core/bitstreams/02df9974-b461-4e1e-b055-70390e7482d9/content

## Citation
If you use this code or method, please cite:

Carrillo-Munoz, Maria, and Bhisham Sharma. “Identifying elastic wave polarization and bandgaps in periodic solid media.” International Journal of Mechanical Sciences 252 (2023): 108363. https://www.sciencedirect.com/science/article/pii/S0020740323002655

## License
MIT License

For questions, contact: Maria Carrillo Munoz
