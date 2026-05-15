"""
plot_apl_contour_definitive.py
==============================
APL contour plot with:
  1. Smooth cubic interpolation over full box
  2. Trajectory-based leaflet assignment (robust)
  3. Protein represented as black CA dots only (Arial Bold Styling)
  4. Red=High APL (Bulk), Blue=Low APL (Interface Compression)
"""

import argparse
import glob
import os
import sys
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Polygon as MplPolygon
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import pandas as pd
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter
from shapely.geometry import MultiPoint

# ── plot parameters ───────────────────────────────────────────────────────────
GRID_N         = 200        
SMOOTH_SIGMA   = 2.5        
CMAP           = "RdBu_r"   # Red=High (Bulk), Blue=Low (Interface)
CONTOUR_LVLS   = 16         
CA_DOT_SIZE    = 6          
CA_DOT_COLOR   = "black"
CA_DOT_ALPHA   = 0.7
DPI            = 200
FIG_WIDTH      = 22         
FIG_HEIGHT     = 7.5

# Global Style Settings for Publication
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Liberation Sans', 'DejaVu Sans']
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'
# ─────────────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--pattern", required=True)
    p.add_argument("--gro",     required=True)
    p.add_argument("--xtc",     required=True)
    p.add_argument("--pdb",     required=True)
    p.add_argument("--out",     required=True)
    p.add_argument("--name",    required=True)
    p.add_argument("--stride",  type=int, default=5)
    return p.parse_args()

def read_gro(gro_path):
    with open(gro_path, "r") as fh:
        lines = fh.readlines()
    box_vals = lines[-1].split()
    box_x = float(box_vals[0])
    box_y = float(box_vals[1])
    ca_xy = []
    for line in lines[2:-1]:
        if len(line) < 44: continue
        if line[10:15].strip() == "CA":
            try: ca_xy.append([float(line[20:28]), float(line[28:36])])
            except: continue
    return box_x, box_y, np.array(ca_xy) if ca_xy else None

def read_pdb_ca(pdb_path):
    ca_xy = []
    if not os.path.exists(pdb_path): return None
    with open(pdb_path, "r") as fh:
        for line in fh:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                try: ca_xy.append([float(line[30:38])/10.0, float(line[38:46])/10.0])
                except: continue
    return np.array(ca_xy) if ca_xy else None

def assign_leaflets_from_traj(gro, xtc, stride):
    import MDAnalysis as mda
    print(f"  Performing trajectory-based leaflet and protein alignment...")
    # Load with bond guessing to allow fragment-aware operations
    u = mda.Universe(gro, xtc, guess_bonds=True)
    
    heads = u.select_atoms("resname POPC and name P")
    prot = u.select_atoms("protein")
    prot_ca = u.select_atoms("protein and name CA")
    
    resids = heads.resids
    unique_resids = np.unique(resids)
    votes = {int(r): 0 for r in unique_resids}
    midplanes = []
    
    # To store centered CA positions across frames
    ca_coords_accum = []
    
    for ts in u.trajectory[::stride]:
        # Center the protein (as an intact fragment) in the box
        prot.wrap(compound='fragments')
        com = prot.center_of_mass()
        u.atoms.translate((ts.dimensions[:3]/2.0) - com)
        u.atoms.wrap(compound='atoms')
        
        # Assignments
        z = heads.positions[:, 2]
        mid = np.mean(z)
        midplanes.append(mid)
        for i, r in enumerate(resids):
            votes[int(r)] += 1 if z[i] > mid else -1
            
        # Collect CA coords (nm)
        ca_coords_accum.append(prot_ca.positions[:, :2] / 10.0)
            
    assignment = {r: ("upper" if v >= 0 else "lower") for r, v in votes.items()}
    # Calculate average CA positions
    avg_ca_coords = np.mean(ca_coords_accum, axis=0)
    
    return assignment, np.mean(midplanes) / 10.0, avg_ca_coords

def make_grid(x, y, apl, box_x, box_y, n=GRID_N):
    lo, hi = np.nanpercentile(apl, [0.5, 99.5])
    mask = (apl >= lo) & (apl <= hi) & np.isfinite(apl)
    x, y, apl = x[mask], y[mask], apl[mask]
    ox, oy = np.array([-box_x, 0, box_x]), np.array([-box_y, 0, box_y])
    xt = np.concatenate([x + dx for dx in ox for dy in oy])
    yt = np.concatenate([y + dy for dx in ox for dy in oy])
    at = np.tile(apl, 9)
    xi = np.linspace(0, box_x, n)
    yi = np.linspace(0, box_y, n)
    Xi, Yi = np.meshgrid(xi, yi)
    grid = griddata(np.column_stack([xt, yt]), at, (Xi, Yi), method="cubic", fill_value=np.nan)
    nan_mask = np.isnan(grid)
    if nan_mask.any():
        grid_nn = griddata(np.column_stack([xt, yt]), at, (Xi, Yi), method="nearest")
        grid[nan_mask] = grid_nn[nan_mask]
    return gaussian_filter(grid, sigma=SMOOTH_SIGMA), xi, yi

def draw_panel(ax, grid, xi, yi, box_x, box_y, ca_xy_plot, title, vmin, vmax):
    X, Y = np.meshgrid(xi, yi)
    cf = ax.contourf(X, Y, grid, levels=np.linspace(vmin, vmax, CONTOUR_LVLS), cmap=CMAP, vmin=vmin, vmax=vmax, extend="both")
    ax.contour(X, Y, grid, levels=np.linspace(vmin, vmax, CONTOUR_LVLS), colors="k", linewidths=0.2, alpha=0.3, vmin=vmin, vmax=vmax)
    
    if ca_xy_plot is not None:
        ax.scatter(ca_xy_plot[:, 0], ca_xy_plot[:, 1], c=CA_DOT_COLOR, s=CA_DOT_SIZE, alpha=CA_DOT_ALPHA, linewidths=0, zorder=5)

    ax.set_xlim(0, box_x); ax.set_ylim(0, box_y)
    ax.set_title(title, fontsize=20, fontweight="bold", pad=10)
    ax.set_xlabel("Box X (nm)", fontsize=18, fontweight="bold")
    ax.set_ylabel("Box Y (nm)", fontsize=18, fontweight="bold")
    ax.tick_params(labelsize=16, width=2.5, length=8)
    for spine in ax.spines.values(): spine.set_linewidth(2.5)
    ax.set_aspect("equal")
    return cf

def process(args):
    print(f"\nProcessing: {args.name}")
    box_x, box_y, _ = read_gro(args.gro)
    
    # Get assignment and PERFECTLY ALIGNED CA COORDS from trajectory
    assignment, midplane_nm, ca_plot = assign_leaflets_from_traj(args.gro, args.xtc, args.stride)
    
    # Note: ca_plot and assignment are calculated using the SAME centering logic
    # as was used in generate_apl_csvs.py (protein COM to box center).
    # So we don't need manual shift sx, sy here as long as we use [0, box] wrapping.

    df = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(args.pattern))[::args.stride]])
    df["leaflet"] = df["resid"].astype(int).map(assignment)
    df = df.dropna(subset=["leaflet"])
    
    # Ensure coordinates are wrapped correctly to the [0, box] range
    df["x_nm"] = df["x_nm"] % box_x
    df["y_nm"] = df["y_nm"] % box_y
    
    ug, xi, yi = make_grid(df[df["leaflet"]=="upper"]["x_nm"].values, df[df["leaflet"]=="upper"]["y_nm"].values, df[df["leaflet"]=="upper"]["apl_nm2"].values, box_x, box_y)
    lg, _, _ = make_grid(df[df["leaflet"]=="lower"]["x_nm"].values, df[df["leaflet"]=="lower"]["y_nm"].values, df[df["leaflet"]=="lower"]["apl_nm2"].values, box_x, box_y)
    avg_g = (ug + lg) / 2.0
    
    all_vals = np.concatenate([ug.ravel(), lg.ravel()])
    all_vals = all_vals[np.isfinite(all_vals)]
    vmin, vmax = np.percentile(all_vals, [2, 98])
    print("  APL range (2nd-98th pct): " + str(round(vmin, 4)) + " - " + str(round(vmax, 4)) + " nm^2")

    fig, axes = plt.subplots(1, 3, figsize=(FIG_WIDTH, FIG_HEIGHT), constrained_layout=True)
    fig.suptitle(f"{args.name} Area Per Lipid", fontsize=32, fontweight="bold", y=1.05)
    
    for ax, grid, title in zip(axes, [ug, lg, avg_g], ["Upper Leaflet", "Lower Leaflet", "Membrane"]):
        cf = draw_panel(ax, grid, xi, yi, box_x, box_y, ca_plot, title, vmin, vmax)
    
    sm = ScalarMappable(cmap=CMAP, norm=Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, orientation="vertical", fraction=0.015, pad=0.01, shrink=0.85)
    cbar.set_label("APL (nm²)", fontsize=20, fontweight="bold")
    cbar.ax.tick_params(labelsize=16)
    
    out_path = os.path.join(args.out, f"{args.name}_Area_Per_Lipid.pdf")
    fig.savefig(out_path, bbox_inches="tight")
    print(f"Saved: {out_path}")

if __name__ == "__main__":
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)
    process(args)
