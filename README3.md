Lattice Dispersion Visualization (MATLAB)
This repository provides a MATLAB script for importing, analyzing, and visualizing the dynamic behavior of 2D rectangular lattice structures from Abaqus eigenfrequency analysis.

Generating the .rpt files (Abaqus or COMSOL)
You can create compatible text reports from either Abaqus (.rpt) or COMSOL (plain-text/CSV exported and saved with the same filenames). Keep it simple and consistent:
• Include the eigenfrequency for each mode at each k-point.  
• Optionally include any modal quantities you have available (e.g., generalized/effective masses, participation factors, rotations). It’s okay if some fields are missing because certain DOFs are constrained—just keep the same set of columns for all k-points in a given file.  
• Use consistent units across all k-points (e.g., frequency in Hz or rad/s; k in 1/m if included).  
• For a PATH file, sweep along your chosen high-symmetry path and concatenate the results in order.  
• For a GRID file, sweep a uniform (kx, ky) grid and concatenate in a consistent ordering (any ordering is fine as long as it’s the same throughout).

Abaqus (concise workflow)
1) Open the .odb in Visualization.  
2) Create XY data from History Output for your eigenfrequency step (Whole Model). Select eigenfrequency and any other modal quantities you wish to export.  
3) Report → XY… to write a text file; append/collect all k-points into one file. Name it `abaqus_Frequencies_PATH.rpt` (path) or `abaqus_Frequencies_GRID.rpt` (grid), or update the script paths accordingly.

Abaqus (scripted hint)
If you automate k-sweeps, use `session.XYDataFromHistory(...)` to gather the series you have (e.g., eigenfrequency plus any available modal terms) and `session.writeXYReport(fileName=..., xyData=(...))` to append them—keep the column order consistent across k-points.

COMSOL (concise workflow)
1) Run an Eigenfrequency study with Bloch/Floquet periodicity; sweep k along a path or over a (kx, ky) grid.  
2) Collect results in a Table (eigenfrequencies and any modal metrics you compute).  
3) Export the Table as space-delimited text or CSV and save with the expected filenames (`abaqus_Frequencies_PATH.rpt` or `abaqus_Frequencies_GRID.rpt`) or adjust the MATLAB script paths. Ensure headers/column order are consistent across all k-points.

Main features:

Path Analysis: Frequency dispersion curves along high-symmetry paths in the Brillouin Zone (using abaqus_Frequencies_PATH.rpt).
Grid Analysis: 3D dispersion surfaces over a grid of wave vectors (kx, ky) (using abaqus_Frequencies_GRID.rpt).
Polarization Visualization: Effective mass data for translation & rotation modes.
Group Velocity Analysis: Direction and speed of energy propagation.
Getting Started
Prepare Abaqus Output Files:
Export eigenfrequency data (including effective mass and rotational mass) as .rpt files.
Required files:

abaqus_Frequencies_PATH.rpt (for boundary/path analysis)
abaqus_Frequencies_GRID.rpt (for full IBZ grid analysis) Place these files in the example_data/ folder.
Edit File Paths (if needed):

In the script, set the paths to your local .rpt files if different.
Run the Script:

Open Rectangular_lattice_3by1.m in MATLAB and run it.
Customize & Extend:

All code is commented for easy extension (e.g., different lattice geometries, advanced plots).
Example Data
Example Abaqus output files are provided in example_data/:

abaqus_Frequencies_PATH.rpt — for path analysis
abaqus_Frequencies_GRID.rpt — for grid analysis
Use these files to test and reproduce the MATLAB analysis.

Note:

The code and .rpt files provided correspond to the rectangular lattice illustrated in Figure 6 of the journal article "Identifying elastic wave polarization and bandgaps in periodic solid media" and to Section 5.3.1 of the dissertation "A Novel Method to Study the Dynamic Behavior of Lattice Structures".

Further reading
For additional background and methodology details, see the dissertation:

Carrillo-Munoz, Maria J. "A Novel Method to Study the Dynamic Behavior of Lattice Structures." PhD Dissertation, Wichita State University, 2023. https://soar.wichita.edu/server/api/core/bitstreams/02df9974-b461-4e1e-b055-70390e7482d9/content
Citation
If you use this code or method, please cite:

Carrillo-Munoz, Maria, and Bhisham Sharma. “Identifying elastic wave polarization and bandgaps in periodic solid media.” International Journal of Mechanical Sciences 252 (2023): 108363. https://www.sciencedirect.com/science/article/pii/S0020740323002655
License
MIT License

For questions, contact: Maria Carrillo Munoz
