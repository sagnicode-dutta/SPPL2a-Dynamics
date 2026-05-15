import MDAnalysis as mda
import numpy as np
import pandas as pd
import os

def calculate_thickness(u):
    """
    Calculates total thickness, upper leaflet, and lower leaflet thickness per frame.
    """
    phosphates = u.select_atoms("resname POPC and name P")
    data = []

    print(f"  Processing {len(u.trajectory[::10])} frames...")
    
    for ts in u.trajectory[::10]:
        com_z = u.select_atoms("resname POPC").center_of_mass()[2]
        
        top_leaflet = phosphates[phosphates.positions[:, 2] > com_z]
        bottom_leaflet = phosphates[phosphates.positions[:, 2] < com_z]
        
        if len(top_leaflet) > 0 and len(bottom_leaflet) > 0:
            z_top = np.mean(top_leaflet.positions[:, 2])
            z_bottom = np.mean(bottom_leaflet.positions[:, 2])
            
            total = z_top - z_bottom
            upper = z_top - com_z
            lower = com_z - z_bottom # Distance is positive
            
            data.append({
                "Time_ns": ts.time / 1000.0,
                "Total": total,
                "Upper": upper,
                "Lower": lower
            })
            
    return pd.DataFrame(data)

def run_all():
    systems = [
        {"name": "9K92", "gro": "9K92_main/step7_1.gro", "xtc": "analysis/9K92/reduced.xtc"},
        {"name": "9K93", "gro": "9K93_main/step7_1.gro", "xtc": "analysis/9K93/reduced.xtc"},
        {"name": "9K95", "gro": "9K95_main/step7_1.gro", "xtc": "analysis/9K95/reduced.xtc"},
    ]

    all_data = []
    for sys_info in systems:
        print(f"Processing {sys_info['name']}...")
        u = mda.Universe(sys_info['gro'], sys_info['xtc'])
        df_sys = calculate_thickness(u)
        df_sys["System"] = sys_info['name']
        all_data.append(df_sys)
    
    final_df = pd.concat(all_data)
    final_df.to_csv("analysis/membrane_thickness_detailed.csv", index=False)
    print("Detailed membrane thickness data saved to analysis/membrane_thickness_detailed.csv")

if __name__ == "__main__":
    run_all()
