import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import numpy as np
import os

# Publication-quality settings
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1.2

# System-specific colors (from project standards)
SYSTEM_COLORS = {
    'Bound SPPL2a': '#53c719',  # Green
    'PS1 (Bound)': '#9929ea'    # Purple
}

# Component colors (distinct palette for the stacks)
COMPONENT_COLORS = {
    'vdW': '#56B4E9',           # Sky blue
    'Elec': '#E69F00',          # Orange  
    'PolarSolv': '#CC79A7',     # Purple
    'NonPolarSolv': '#009E73'   # Green
}

def load_decomp_data(filepath, system_name, exclude_list=None):
    """Parses the MM/PBSA decomposition file and filters results."""
    if exclude_list is None:
        exclude_list = []
        
    data = []
    if not os.path.exists(filepath):
        print(f"Warning: {filepath} not found.")
        return pd.DataFrame()

    with open(filepath, 'r') as f:
        in_data = False
        for line in f:
            if 'Total Energy Decomposition:' in line:
                in_data = True
                continue
            # Exclude lines starting with L: (Ligand) and use explicit exclusion list
            if in_data and line.startswith('R:'):
                parts = line.strip().split(',')
                if len(parts) >= 18:
                    res_info = parts[0].split(':')
                    resname = res_info[2]
                    resid = res_info[3]
                    label = f"{resname}{resid}"
                    
                    if label in exclude_list:
                        continue
                        
                    data.append({
                        'Residue': label,
                        'vdW': float(parts[4]),
                        'Elec': float(parts[7]),
                        'PolarSolv': float(parts[10]),
                        'NonPolarSolv': float(parts[13]),
                        'Total': float(parts[16]),
                        'TotalSEM': float(parts[18]),
                        'System': system_name
                    })
            if in_data and 'Sidechain Energy Decomposition:' in line:
                break
                
    df = pd.DataFrame(data)
    if df.empty:
        return df
        
    # Filter by threshold (-5 kcal/mol) and sort
    df = df[df['Total'] < -5.0].sort_values('Total', ascending=True)
    return df

def plot_stacked_components(ax, df, title, system_color):
    """Creates a horizontal stacked bar plot for a single system."""
    if df.empty:
        ax.text(0.5, 0.5, "No data available", ha='center')
        return

    # Reverse for plotting (best contributors at top)
    df_plot = df.iloc[::-1].head(12) # Show top 12 contributors
    y_pos = np.arange(len(df_plot))
    
    # Plot favorable (negative) components
    cumulative_left = np.zeros(len(df_plot))
    for comp in ['vdW', 'Elec', 'NonPolarSolv']:
        vals = np.minimum(df_plot[comp].values, 0)
        ax.barh(y_pos, vals, left=cumulative_left, color=COMPONENT_COLORS[comp], 
                edgecolor='white', linewidth=0.5)
        cumulative_left += vals
        
    # Plot unfavorable (positive) components (Polar Solvation usually)
    vals_pos = np.maximum(df_plot['PolarSolv'].values, 0)
    ax.barh(y_pos, vals_pos, color=COMPONENT_COLORS['PolarSolv'], 
            edgecolor='white', linewidth=0.5)
    
    # Add "Total Energy" dots with the specific system color
    ax.plot(df_plot['Total'], y_pos, 'o', color=system_color, markersize=8, 
            markeredgecolor='black', markeredgewidth=1, label='Total ΔG')
    
    # Customization
    ax.set_yticks(y_pos)
    ax.set_yticklabels(df_plot['Residue'], fontweight='bold')
    ax.set_title(title, fontweight='bold', color=system_color, fontsize=12, loc='left')
    ax.axvline(x=0, color='black', linewidth=1.2, alpha=0.5)
    ax.grid(axis='x', alpha=0.2, linestyle='--')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel('Contribution (kcal/mol)')

# Load data with exclusions
# PS1: Exclude Asp169 artifact (if found), and highly flexible termini
# Note: load_decomp_data now automatically skips L: (Ligands)
ps1_df = load_decomp_data("analysis/mmpbsa/ps1_mmpbsa/FINAL_DECOMP_MMPBSA_PS1_v2.dat", 
                         "PS1 (Bound)", exclude_list=['ASP169', 'ARG269', 'GLU273'])
sppl2a_df = load_decomp_data("analysis/mmpbsa/sppl2a_mmpbsa/FINAL_DECOMP_MMPBSA_SPPL2A_v2.dat", 
                            "Bound SPPL2a")

# Create figure
fig = plt.figure(figsize=(12, 8))
gs = GridSpec(1, 2, figure=fig, wspace=0.35)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])

# Plot both systems
plot_stacked_components(ax1, ps1_df, "A) PS1 (Bound)", SYSTEM_COLORS['PS1 (Bound)'])
plot_stacked_components(ax2, sppl2a_df, "B) Bound SPPL2a", SYSTEM_COLORS['Bound SPPL2a'])

# Add shared legend
legend_handles = [
    mpatches.Patch(color=COMPONENT_COLORS['vdW'], label='vdW'),
    mpatches.Patch(color=COMPONENT_COLORS['Elec'], label='Electrostatic'),
    mpatches.Patch(color=COMPONENT_COLORS['NonPolarSolv'], label='Nonpolar Solv.'),
    mpatches.Patch(color=COMPONENT_COLORS['PolarSolv'], label='Polar Solv.'),
    plt.Line2D([0], [0], marker='o', color='w', label='Total Binding Energy',
               markerfacecolor='gray', markersize=10, markeredgecolor='black')
]
fig.legend(handles=legend_handles, loc='upper center', bbox_to_anchor=(0.5, 1.02),
           ncol=5, frameon=False, fontsize=10)

plt.tight_layout()
# Create plots directory if it doesn't exist
os.makedirs("analysis/plots", exist_ok=True)
plt.savefig("analysis/plots/mmpbsa_decomposition_publication.pdf", bbox_inches='tight', dpi=300)
plt.savefig("analysis/plots/mmpbsa_decomposition_publication.png", bbox_inches='tight', dpi=300)

print("Publication-quality plots generated: analysis/plots/mmpbsa_decomposition_publication.pdf")
