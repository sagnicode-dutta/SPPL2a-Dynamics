"""
generate_apl_csvs.py
====================
Generates per-frame APL CSVs using protein-aware Voronoi tessellation.

The fix vs the original script:
  - Protein CA atoms are inserted as ghost seeds into the Voronoi calculation.
  - Their cells absorb the protein-occupied area.
  - Ghost cell areas are discarded — not assigned to any lipid.
  - Result: lipid cells at the protein interface are correctly bounded,
    no artificial expansion into the protein core, no blue compression ring.
"""

import MDAnalysis as mda
import numpy as np
import pandas as pd
from scipy.spatial import Voronoi
from shapely.geometry import Polygon, box as sbox
import os
import warnings
warnings.filterwarnings("ignore")


# ── CONFIG ────────────────────────────────────────────────────────────────────

SYSTEMS = [
    {"name": "9K92", "gro": "9K92_main/step7_1.gro", "xtc": "analysis/9K92/reduced.xtc"},
    {"name": "9K93", "gro": "9K93_main/step7_1.gro", "xtc": "analysis/9K93/reduced.xtc"},
    {"name": "9K95", "gro": "9K95_main/step7_1.gro", "xtc": "analysis/9K95/reduced.xtc"},
]

STRIDE    = 10
OUT_BASE  = "analysis/apl_raw"
LIPID     = "POPC"
HEAD_ATOM = "P"

# Max plausible APL (nm^2). Cells larger than this are clipped to NaN.
# Typical POPC APL ~ 0.65 nm^2. 3x that is a safe ceiling.
APL_MAX = 2.0
APL_MIN = 0.2

# ─────────────────────────────────────────────────────────────────────────────


def voronoi_apl_protein_aware(lipid_xy, protein_xy, box_dims):
    """
    Compute Voronoi APL for each lipid headgroup, correctly bounded by protein.

    Parameters
    ----------
    lipid_xy   : (N, 2) array, lipid headgroup XY in nm, wrapped to [0,Lx]x[0,Ly]
    protein_xy : (M, 2) array, protein CA XY in nm (ghost seeds), or None
    box_dims   : (Lx, Ly) in nm

    Returns
    -------
    areas : (N,) array of APL in nm^2. NaN for degenerate cells.
    """
    Lx, Ly = box_dims

    # Combine lipid + protein ghost seeds
    if protein_xy is not None and len(protein_xy) > 0:
        all_xy   = np.vstack([lipid_xy, protein_xy])
        n_lipid  = len(lipid_xy)
        n_ghost  = len(protein_xy)
    else:
        all_xy   = lipid_xy
        n_lipid  = len(lipid_xy)
        n_ghost  = 0

    N_total = len(all_xy)

    # 3x3 PBC replication
    offsets = np.array([
        [dx * Lx, dy * Ly]
        for dx in [-1, 0, 1]
        for dy in [-1, 0, 1]
    ])
    tiled = np.vstack([all_xy + off for off in offsets])

    try:
        vor = Voronoi(tiled)
    except Exception:
        return np.full(n_lipid, np.nan)

    central_start = 4 * N_total   # index of first central-replica point
    box_poly = sbox(0, 0, Lx, Ly)

    areas = np.full(n_lipid, np.nan)

    for i in range(n_lipid):
        ridx   = vor.point_region[central_start + i]
        region = vor.regions[ridx]

        if -1 in region or len(region) == 0:
            continue

        verts = vor.vertices[region]
        try:
            cell = Polygon(verts)
        except Exception:
            continue

        # Clip to simulation box
        cell = cell.intersection(box_poly)
        if cell.is_empty:
            continue

        areas[i] = cell.area   # nm^2

    return areas


def run_system(name, gro, xtc):
    print("\n" + "=" * 50)
    print("Processing: " + name)
    print("=" * 50)

    u = mda.Universe(gro, xtc)

    heads = u.select_atoms("resname " + LIPID + " and name " + HEAD_ATOM)
    if len(heads) == 0:
        print("ERROR: No headgroup atoms found. Check LIPID/HEAD_ATOM.")
        return

    try:
        protein_ca = u.select_atoms("protein and name CA")
        if len(protein_ca) == 0:
            protein_ca = u.select_atoms("protein")
        has_protein = len(protein_ca) > 0
    except Exception:
        has_protein = False

    print("  Lipid P atoms : " + str(len(heads)))
    print("  Protein CA    : " + str(len(protein_ca) if has_protein else 0))

    out_dir = os.path.join(OUT_BASE, name)
    os.makedirs(out_dir, exist_ok=True)

    traj = list(u.trajectory[::STRIDE])
    n_frames = len(traj)
    print("  Frames        : " + str(n_frames))

    for fi, ts in enumerate(traj):
        if fi % max(1, n_frames // 10) == 0:
            pct = int(100 * fi / n_frames)
            print("  Frame " + str(fi + 1) + "/" + str(n_frames) + " (" + str(pct) + "%)")

        box    = ts.dimensions[:3] / 10.0   # nm
        Lx, Ly = box[0], box[1]

        # Lipid headgroup positions, PBC-wrapped
        lip_pos = heads.positions / 10.0
        lip_pos = lip_pos % np.array([Lx, Ly, box[2]])

        # Protein CA positions, PBC-wrapped
        if has_protein:
            prot_pos = protein_ca.positions / 10.0
            prot_xy  = prot_pos[:, :2] % np.array([Lx, Ly])
        else:
            prot_xy = None

        # Leaflet assignment: z relative to membrane midplane
        z_mid       = float(np.median(lip_pos[:, 2]))
        upper_mask  = lip_pos[:, 2] > z_mid
        lower_mask  = ~upper_mask

        frame_data = []

        for leaf_name, mask in [("upper", upper_mask), ("lower", lower_mask)]:
            sel_pos    = lip_pos[mask]
            sel_resids = heads.resids[mask]

            if len(sel_pos) < 4:
                continue

            sel_xy = sel_pos[:, :2]

            # ── PROTEIN-AWARE VORONOI ─────────────────────────────────────
            areas = voronoi_apl_protein_aware(sel_xy, prot_xy, (Lx, Ly))
            # ─────────────────────────────────────────────────────────────

            for j in range(len(sel_pos)):
                apl = areas[j]
                # Sanity-clip: discard physically impossible values
                if np.isfinite(apl) and (apl < APL_MIN or apl > APL_MAX):
                    apl = np.nan
                frame_data.append({
                    "resid":   int(sel_resids[j]),
                    "leaflet": leaf_name,
                    "x_nm":    round(float(sel_pos[j, 0]), 4),
                    "y_nm":    round(float(sel_pos[j, 1]), 4),
                    "apl_nm2": round(float(apl), 6) if np.isfinite(apl) else np.nan,
                })

        df = pd.DataFrame(frame_data)
        df.to_csv(
            os.path.join(out_dir, "apl_raw_frame_" + str(fi) + ".csv"),
            index=False
        )

    print("  Done: " + out_dir)


if __name__ == "__main__":
    for sys_cfg in SYSTEMS:
        try:
            run_system(sys_cfg["name"], sys_cfg["gro"], sys_cfg["xtc"])
        except Exception as exc:
            import traceback
            print("ERROR in " + sys_cfg["name"] + ": " + str(exc))
            traceback.print_exc()

    print("\nAll systems complete.")
