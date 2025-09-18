# Lattice Dispersion Visualization (MATLAB)
This repository provides a MATLAB script for importing, analyzing, and visualizing the dynamic behavior of 2D rectangular lattice structures from Abaqus eigenfrequency analysis.

## Generating the `.rpt` files (Abaqus & COMSOL)
You can produce compatible text reports either directly from **Abaqus** (`.rpt`) or from **COMSOL** (plain-text export with the same column names/order). The MATLAB script expects, per k-point, the following series (in this order), which appear in the sample `abaqus_Frequencies_*.rpt` files:

1) `Eigenfrequency: EIGFREQ for Whole Model`  
2) `Generalized mass: GM for Whole Model`  
3) `Effective mass, x-component`  
4) `Effective mass, y-component`  
5) `Effective mass, z-component`  
6) `Effective mass, x-rotation`  
7) `Effective mass, y-rotation`  
8) `Effective mass, z-rotation`  
9) `Participation factor, x-component`  
10) `Participation factor, y-component`  
11) `Participation factor, z-component`  
12) `Participation factor, x-rotation`  
13) `Participation factor, y-rotation`  
14) `Participation factor, z-rotation`

> **PATH file** (`abaqus_Frequencies_PATH.rpt`): concatenate results along your chosen high-symmetry path.  
> **GRID file** (`abaqus_Frequencies_GRID.rpt`): concatenate results over a uniform `(kx, ky)` grid (row-major order is fine, but be consistent).

### Abaqus (GUI, concise)
1. **Open your `.odb`** in *Visualization*.  
2. **Create XY data from History Output** for the eigenfrequency step of each k-point (Whole Model):  
   - *Eigenfrequency (EIGFREQ)* and *Generalized mass (GM)*.  
   - *Effective mass* (x, y, z **components** and x, y, z **rotations**).  
   - *Participation factors* (PF1..PF3 components and PF4..PF6 rotations).  
3. **Report to file**: *Report → XY…* and select the XY objects **in the order listed above**.  
   - Save as `example_data/abaqus_Frequencies_PATH.rpt` (for path) or `example_data/abaqus_Frequencies_GRID.rpt` (for grid).

### Abaqus (scripted, brief)
If you sweep k-points programmatically, you can automate report generation with:
- `session.XYDataFromHistory(...)` to collect each series for a given step/k-point, and  
- `session.writeXYReport(fileName='abaqus_Frequencies_GRID.rpt', xyData=(...))`  
exactly as in your Python snippet (the order of `xyData=(...)` should match the list above).

### COMSOL (plain-text export)
1. **Eigenfrequency study** on the unit cell with **Bloch/Floquet periodicity**; sweep k along a **path** (single parameter) or over a **grid** (two parameters).  
2. Create a **Table** that includes the **eigenfrequency** and any **global evaluations** you use for generalized/effective masses and participation factors (if you compute them).  
3. **Export → Data** the Table as a **space-delimited text** (or CSV) file.  
4. **Rename headers** to match the list above (same wording) and save as:  
   - `example_data/abaqus_Frequencies_PATH.rpt` (path) or  
   - `example_data/abaqus_Frequencies_GRID.rpt` (grid).  
   Keep units consistent across all k-points.

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
