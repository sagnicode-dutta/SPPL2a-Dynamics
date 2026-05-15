import MDAnalysis as mda
import numpy as np
import pandas as pd
import os
import sys

def calculate_scd(u, chain_prefix, carbon_range, h_suffix_map):
    scd_data = {i: [0.0, 0.0] for i in carbon_range}
    
    # Pre-identify atom groups for each carbon and hydrogen type
    # This avoids repeated selections
    c_groups = {}
    h_groups = {}
    
    for i in carbon_range:
        c_name = f"{chain_prefix}{i}"
        c_groups[i] = u.select_atoms(f"resname POPC and name {c_name}")
        
        h_suffixes = h_suffix_map.get(i, [])
        h_groups[i] = {}
        for h_suffix in h_suffixes:
            h_name = f"H{i}{h_suffix}"
            h_groups[i][h_suffix] = u.select_atoms(f"resname POPC and name {h_name}")

    print(f"  Processing {len(u.trajectory[::10])} frames...")
    
    for ts in u.trajectory[::10]:
        z_axis = np.array([0, 0, 1])
        for i in carbon_range:
            c_atoms = c_groups[i]
            for h_suffix, h_atoms in h_groups[i].items():
                if len(c_atoms) == len(h_atoms) and len(c_atoms) > 0:
                    # Vectorized calculation
                    # Ensure they are aligned by residue
                    # (Assuming standard ordering in the topology)
                    vecs = h_atoms.positions - c_atoms.positions
                    norms = np.linalg.norm(vecs, axis=1)
                    cos_thetas = vecs[:, 2] / norms # Dot product with [0,0,1]
                    scds = -0.5 * (3 * cos_thetas**2 - 1)
                    
                    scd_data[i][0] += np.sum(scds)
                    scd_data[i][1] += len(scds)
    
    final_scd = []
    for i in carbon_range:
        if scd_data[i][1] > 0:
            avg_scd = scd_data[i][0] / scd_data[i][1]
            # Scd is traditionally reported as positive |Scd|
            final_scd.append({"Carbon": i, "SCD": abs(avg_scd)})
    
    return pd.DataFrame(final_scd)
                            
    # Final averaging
    final_scd = []
    for i in carbon_range:
        if scd_data[i][1] > 0:
            avg_scd = scd_data[i][0] / scd_data[i][1]
            final_scd.append({"Carbon": i, "SCD": avg_scd})
    
    return pd.DataFrame(final_scd)

def run_all():
    systems = [
        {"name": "9K92", "gro": "9K92_main/step7_1.gro", "xtc": "analysis/9K92/reduced.xtc"},
        {"name": "9K93", "gro": "9K93_main/step7_1.gro", "xtc": "analysis/9K93/reduced.xtc"},
        {"name": "9K95", "gro": "9K95_main/step7_1.gro", "xtc": "analysis/9K95/reduced.xtc"},
    ]

    # Map for sn-2 (Oleoyl, C2)
    # C22-C28: R, S
    # C29, C210: 1
    # C211-C218: R, S
    c2_h_map = {i: ['R', 'S'] for i in range(2, 19)}
    c2_h_map[9] = ['1']
    c2_h_map[10] = ['1']

    # Map for sn-1 (Palmitoyl, C3)
    # C32-C316: X, Y
    c3_h_map = {i: ['X', 'Y'] for i in range(2, 17)}

    for sys_info in systems:
        print(f"Processing {sys_info['name']}...")
        u = mda.Universe(sys_info['gro'], sys_info['xtc'])
        
        print("  Calculating sn-1 (Palmitoyl)...")
        df_sn1 = calculate_scd(u, "C3", range(2, 17), c3_h_map)
        df_sn1.to_csv(f"analysis/{sys_info['name']}_scd_palmitoyl.csv", index=False)
        
        print("  Calculating sn-2 (Oleoyl)...")
        df_sn2 = calculate_scd(u, "C2", range(2, 19), c2_h_map)
        df_sn2.to_csv(f"analysis/{sys_info['name']}_scd_oleoyl.csv", index=False)

if __name__ == "__main__":
    run_all()
