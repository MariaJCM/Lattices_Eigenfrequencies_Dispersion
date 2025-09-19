# Lattice Dispersion Visualization (MATLAB)
This repository provides a MATLAB script for importing, analyzing, and visualizing the dynamic behavior of 2D rectangular lattice structures from Abaqus eigenfrequency analysis.

## Generating input files (.rpt /.txt /.csv) from Abaqus or COMSOL

You can create the required input as **plain-text tables** from either **Abaqus** (e.g., `.rpt`, `.txt`) or **COMSOL** (e.g., `.txt`, `.csv`). The script only needs your **eigenfrequencies** and, the **effective mass components** (translation and/or rotation) for the **DOFs active in your model**. If your export uses different header names or column order, simply identify/match them in the MATLAB script.

### Minimum content & formatting
- Eigenfrequency results **per mode** at each **k-point**.
- A consistent **k-point ordering**:  
  - **PATH** file: concatenate k-points along the chosen high-symmetry path.  
  - **GRID** file: k-points laid out on a uniform `(kx, ky)` grid (if you want to study the complete Irreduzible Brillouin Zone).
- Generalized/effective masses, participation factors.

### Abaqus workflow
1. Run the **eigenfrequency** step(s) for each k-point (e.g., Bloch/Floquet periodic conditions).
2. In **Visualization**, create **XY Data from History Output** for the k-point step and select the results you wish to export (at minimum, eigenfrequencies; optionally other quantities if present).
3. Use **Report → XY…** to write a **single text report** per case (commonly `.rpt`; plain text).  
   Save as `Frequencies_PATH.rpt` (path sweep) or `Frequencies_GRID.rpt` (grid sweep).  
   Keep the **same column order** across all k-points.  
   *(Scripting via the Abaqus Python API can also write directly to txt/CSV if you prefer.)*

### COMSOL workflow
1. Set up an **Eigenfrequency** study with **Floquet periodicity** and sweep a **path** in k-space or a **(kx, ky)** grid.
2. Create a **Table** containing eigenfrequencies and the measures you want to include.
3. **Export → Data** the table as **space-delimited text** or **CSV**, keep headers, and save (or rename) as:  
   `Frequencies_PATH.rpt` (path) or `Frequencies_GRID.rpt` (grid).

---

## Main features:
- **Path Analysis:** Frequency dispersion curves along high-symmetry paths in the Brillouin Zone (using `abaqus_Frequencies_PATH.rpt`).
- **Grid Analysis:** Frequency dispersion surfaces over a grid of wave vectors (kx, ky) (using `abaqus_Frequencies_GRID.rpt`).
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

- Carrillo-Munoz, Maria J. "A Novel Method to Study the Dynamic Behavior of Lattice Structures." PhD Dissertation, Wichita State University, 2023. https://soar.wichita.edu/server/api/core/bitstreams/02df9974-b461-4e1e-b055-70390e7482d9/content

## Citation
If you use this code or method, please cite:

- Carrillo-Munoz, Maria, and Bhisham Sharma. “Identifying elastic wave polarization and bandgaps in periodic solid media.” International Journal of Mechanical Sciences 252 (2023): 108363. https://www.sciencedirect.com/science/article/pii/S0020740323002655

## License
MIT License

For questions, contact: Maria Carrillo Munoz
