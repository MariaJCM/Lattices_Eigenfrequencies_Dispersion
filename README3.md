# Lattice Dispersion Visualization (MATLAB)
This repository provides a MATLAB script for importing, analyzing, and visualizing the dynamic behavior of 2D rectangular lattice structures from Abaqus eigenfrequency analysis.

## Generating the `.rpt`/text files (Abaqus & COMSOL)
You may create the required input files either directly from **Abaqus** (`.rpt`) or export **plain-text/CSV** from **COMSOL** and save/rename it as `.rpt`. The key is consistency across all k-points and files. Depending on boundary conditions (e.g., constrained DOFs), some quantities may be absent—include whatever is available, but keep headers, units, and ordering consistent.

### Minimum content & formatting
- Eigenfrequency results **per mode** at each **k-point**.  
- A consistent **k-point ordering**:  
  - **PATH** file: concatenate k-points along the chosen high-symmetry path.  
  - **GRID** file: k-points laid out on a uniform `(kx, ky)` grid (any consistent row/column order is fine).  
- Plain-text with a stable delimiter (space/CSV) and a repeatable header.  
- Units kept consistent across files (e.g., frequency in Hz; wave-vector components in 1/m if included).  
- Optional additional dynamic quantities (e.g., generalized/effective masses, participation factors) may be included when available.

### Abaqus (concise workflow)
1. Run the **eigenfrequency** step(s) for each k-point (e.g., Bloch/Floquet periodic conditions on the unit cell).
2. In **Visualization**, create **XY Data from History Output** for the k-point step and select the results you wish to export (at minimum, eigenfrequencies; optionally other dynamic quantities if present).
3. Use **Report → XY…** to write a **single** `.rpt` file per case:  
   - `abaqus_Frequencies_PATH.rpt` for the path sweep.  
   - `abaqus_Frequencies_GRID.rpt` for the grid sweep.  
   Ensure the same column order across all k-points.

### COMSOL (concise workflow)
1. Set up an **Eigenfrequency** study on the unit cell with **Floquet periodicity** and sweep either a **path** in k-space or a **(kx, ky)** grid.
2. Create a **Table** containing eigenfrequencies (and any other available/global measures you wish to include).
3. **Export → Data** the Table as **space-delimited text** or **CSV**, keep consistent headers/units, and save (or rename) as:  
   - `abaqus_Frequencies_PATH.rpt` for the path sweep.  
   - `abaqus_Frequencies_GRID.rpt` for the grid sweep.

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

- Carrillo-Munoz, Maria J. "A Novel Method to Study the Dynamic Behavior of Lattice Structures." PhD Dissertation, Wichita State University, 2023. https://soar.wichita.edu/server/api/core/bitstreams/02df9974-b461-4e1e-b055-70390e7482d9/content

## Citation
If you use this code or method, please cite:

- Carrillo-Munoz, Maria, and Bhisham Sharma. “Identifying elastic wave polarization and bandgaps in periodic solid media.” International Journal of Mechanical Sciences 252 (2023): 108363. https://www.sciencedirect.com/science/article/pii/S0020740323002655

## License
MIT License

For questions, contact: Maria Carrillo Munoz
