#!/usr/bin/env python3
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

# Mapping of internal resnames to user-friendly names
LIGAND_NAMES = {
    "FTO": "L685,458",
    "A1D6": "Compound E"
}

def analyze_hbonds(gro, xtc, system_name, internal_resname, user_ligand_name, output_csv):
    print(f"Analyzing {system_name} (Ligand: {user_ligand_name})...")
    u = mda.Universe(gro, xtc)
    
    valid_elements = ('O', 'N', 'F', 'S')
    hydrogens_sel = f"(protein or resname {internal_resname}) and (name H* or name h*)"

    h = HydrogenBondAnalysis(
        u,
        donors_sel=f"protein or resname {internal_resname}",
        acceptors_sel=f"protein or resname {internal_resname}",
        hydrogens_sel=hydrogens_sel,
        d_h_a_angle_cutoff=150,
        d_a_cutoff=3.5
    )
    
    h.run(step=10)
    hb_array = h.results.hbonds
    
    pair_frames = {}
    
    for rec in hb_array:
        frame = int(rec[0])
        d_idx = int(rec[1])
        a_idx = int(rec[3])
        
        d_atom = u.atoms[d_idx]
        a_atom = u.atoms[a_idx]
        
        d_is_ligand = d_atom.resname == internal_resname
        a_is_ligand = a_atom.resname == internal_resname
        
        if (d_is_ligand != a_is_ligand):
            d_type = d_atom.name[0]
            a_type = a_atom.name[0]
            
            if d_type in valid_elements and a_type in valid_elements:
                # Donor label (Residue Name Only)
                d_label = user_ligand_name if d_is_ligand else f"{d_atom.resname}{d_atom.resid}"
                
                # Acceptor label (Residue Name Only)
                a_label = user_ligand_name if a_is_ligand else f"{a_atom.resname}{a_atom.resid}"
                
                key = (d_label, a_label)
                if key not in pair_frames:
                    pair_frames[key] = set()
                pair_frames[key].add(frame)
            
    n_analyzed = len(u.trajectory[::10])
    results = []
    for (d_label, a_label), frameset in pair_frames.items():
        occ = (len(frameset) / n_analyzed) * 100
        results.append({
            "Donor": d_label,
            "Acceptor": a_label,
            "Occupancy": occ,
            "System": system_name
        })
        
    df = pd.DataFrame(results)
    if not df.empty:
        df.to_csv(output_csv, index=False)
    return df

def plot_bubble(df, system_name, system_color, output_pdf):
    # Sort and take top interactions
    top_data = df.sort_values(by="Occupancy", ascending=False).head(20)
    
    # Font settings
    font_settings = {'family': 'sans-serif', 'weight': 'bold', 'size': 14}
    plt.rc('font', **font_settings)

    plt.figure(figsize=(12, 10))
    
    # Create custom colormap from white to the system color
    # This ensures the "correct" color gradient for each system
    cmap = LinearSegmentedColormap.from_list("custom_cmap", ["#ffffff", system_color])
    
    # Bubble plot
    scatter = plt.scatter(top_data["Donor"], top_data["Acceptor"], 
                         s=top_data["Occupancy"]*40, # Increased scale for better visibility
                         c=top_data["Occupancy"], 
                         cmap=cmap,
                         alpha=0.8, edgecolors='black', linewidth=1.5)
    
    plt.xlabel("Donor Residue", fontsize=18, fontweight='bold', labelpad=15)
    plt.ylabel("Acceptor Residue", fontsize=18, fontweight='bold', labelpad=15)
    plt.title(f"Protein-Ligand H-Bond Occupancy Bubble Plot\n{system_name}", fontsize=20, fontweight='bold', pad=20)
    
    plt.xticks(rotation=45, ha='right', fontsize=12, fontweight='bold')
    plt.yticks(fontsize=12, fontweight='bold')
    
    # Colorbar
    cbar = plt.colorbar(scatter)
    cbar.set_label("Occupancy (%)", fontsize=16, fontweight='bold')
    cbar.ax.tick_params(labelsize=14)

    # Size Legend - Highly spaced out
    handles, labels = [], []
    for val in [20, 40, 60, 80, 100]:
        handles.append(plt.scatter([], [], s=val*40, c=system_color, alpha=0.6, edgecolors='black'))
        labels.append(f"{val}%")
    
    legend = plt.legend(handles, labels, title="Occupancy (%)", loc="upper left", 
                        bbox_to_anchor=(1.25, 1),
                        fontsize=12, title_fontsize=14, frameon=True,
                        labelspacing=3.0, handletextpad=2.0)
    plt.setp(legend.get_title(), fontweight='bold')

    plt.tight_layout()
    plt.savefig(output_pdf, bbox_inches='tight')
    plt.close()

# Systems and their established colors
# Apo SPPL2a: #ff0046 | Bound SPPL2a (9K93): #53c719 | Bound PS1 (9K95): #9929ea
systems_data = [
    ("9K93_main/step7_1.gro", "analysis/9K93/reduced.xtc", "Bound SPPL2a", "#53c719", "FTO", "L685,458", "analysis/9K93/hbond_ligand_occupancy.csv", "analysis/plots/hbond_ligand_bubble_9K93.pdf"),
    ("9K95_main/step7_1.gro", "analysis/9K95/reduced.xtc", "Bound PS1", "#9929ea", "A1D6", "Compound E", "analysis/9K95/hbond_ligand_occupancy.csv", "analysis/plots/hbond_ligand_bubble_9K95.pdf")
]

os.makedirs("analysis/plots", exist_ok=True)

for gro, xtc, name, color, int_res, user_res, csv_out, pdf_out in systems_data:
    if os.path.exists(gro) and os.path.exists(xtc):
        df = analyze_hbonds(gro, xtc, name, int_res, user_res, csv_out)
        if not df.empty:
            plot_bubble(df, name, color, pdf_out)
    else:
        print(f"Skipping {name}: Files not found.")

print("Protein-Ligand H-bond occupancy bubble plot analysis complete.")