# Rectangular Lattice Dispersion Visualization (MATLAB)

This repository provides a MATLAB script for importing, analyzing, and visualizing the dynamic behavior of 2D rectangular lattice structures from Abaqus eigenfrequency analysis.

**Main features:**
- **Path Analysis:** Frequency dispersion curves along high-symmetry paths in the Brillouin Zone (using `abaqus_Frequencies_PATH.rpt`).
- **Grid Analysis:** 3D dispersion surfaces over a grid of wave vectors (kx, ky) (using `abaqus_Frequencies_GRID.rpt`).
- **Polarization Visualization:** Effective mass data for translation & rotation modes.
- **Group Velocity Analysis:** Direction and speed of energy propagation.

## Getting Started

1. **Prepare Abaqus Output Files:**  
   Export eigenfrequency data (including effective mass and rotational mass) as `.rpt` files.  
   Required files:
   - `abaqus_Frequencies_PATH.rpt` (for boundary/path analysis)
   - `abaqus_Frequencies_GRID.rpt` (for full IBZ grid analysis)
   Place these files in the `example_data/` folder.

2. **Edit File Paths (if needed):**  
   - In the script, set the paths to your local `.rpt` files if different.

3. **Run the Script:**  
   - Open `Rectangular_lattice_3by1.m` in MATLAB and run it.

4. **Customize & Extend:**  
   - All code is commented for easy extension (e.g., different lattice geometries, advanced plots).

## Example Data

Example Abaqus output files are provided in [`example_data/`](example_data):
- [`abaqus_Frequencies_PATH.rpt`](example_data/abaqus_Frequencies_PATH.rpt) — for path analysis
- [`abaqus_Frequencies_GRID.rpt`](example_data/abaqus_Frequencies_GRID.rpt) — for grid analysis

Use these files to test and reproduce the MATLAB analysis.

## Citation

If you use this code or method, please cite:

- Carrillo-Munoz, Maria J. "A Novel Method to Study the Dynamic Behavior of Lattice Structures." PhD diss., Wichita State University, 2023.
- Carrillo-Munoz, Maria, and Bhisham Sharma. "Identifying elastic wave polarization and bandgaps in periodic solid media." *International Journal of Mechanical Sciences* 252 (2023): 108363.

## License

MIT License

---

For questions, contact: Maria Carrillo Munoz
