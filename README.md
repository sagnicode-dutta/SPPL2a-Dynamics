# SPPL2a Dynamics Analysis

This repository contains the full analysis and plotting pipeline for the molecular dynamics study of **SPPL2a** in three states: 
1. **Apo** (PDB: 9K92)
2. **Inhibitor-bound** (PDB: 9K93)
3. **Substrate-bound** (PDB: 9K95 - $\gamma$-secretase complex)

## Repository Structure

The scripts are organized into a logical workflow:

### `01_trajectory_prep/`
Scripts to clean and prepare raw simulation data (removing PBC, centering, and fitting).
* `step0_prep_trajectories.sh`: Generates protein-only trajectories.

### `02_global_metrics/`
Core calculations for protein stability and equilibration.
* `run_analysis.sh`: GROMACS commands for RMSD, Radius of Gyration (Rg), and SASA.
* `calculate_convergence_points.py`: Statistical determination of simulation equilibration.
* `plot_sppl2a_global_metrics.py`: Combined plotting for stability metrics.

### `03_lipid_analysis/`
Analysis of the protein-membrane environment.
* `calculate_membrane_thickness.py` / `plot_membrane_thickness.py`
* `calculate_scd_mda.py` / `plot_lipid_order.py`: Lipid tail order parameters.
* `generate_apl_csvs.py` / `plot_apl_contour_definitive.py`: Area Per Lipid (APL) maps.

### `04_protein_analysis/`
Detailed protein dynamics and interactions.
* **DCCM:** `plot_dccm.R` (Correlated motions).
* **DSSP:** `tidy_dssp.py` and `plot_dssp_heatmap.R` (Secondary structure evolution).
* **Contacts:** `tm2_contact_occupancy.py` and `plot_tm_contacts.py` (Residue contact analysis).
* **H-bonds:** `hbond_occupancy_analysis.py` and `plot_hbond_final.py` (Hydrogen bond networks).

### `05_structural_selection/`
Methods for identifying representative states and extracting key frames.
* `robust_cluster.sh`: RMSD-based clustering to find the dominant structural centroid.
* `perform_snapshot_extraction.sh`: Automated extraction of time-point snapshots (e.g., every 200ns).

### `06_mmpbsa/`
Binding energetics and energy decomposition.
* `run_mmpbsa.sh`: Template script for `gmx_MMPBSA` calculations and input generation.
* `plot_mmpbsa_publication.py`: Main binding energy comparison plots.
* `plot_mmpbsa_decomposition.R`: Residue-wise energy decomposition (vdW vs. Electrostatic).

## Data Availability
The reduced trajectories and structure files required to run these scripts are available on Zenodo: `[Insert Zenodo DOI here]`

## Requirements
- **GROMACS** (Tested on 2023+)
- **gmx_MMPBSA** (For folder 06)
- **Python 3.x** (Required libraries: `MDAnalysis`, `NumPy`, `Pandas`, `Matplotlib`, `Seaborn`)
- **R** (Required libraries: `ggplot2`, `dplyr`, `bio3d`)
