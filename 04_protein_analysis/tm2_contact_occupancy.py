import MDAnalysis as mda
from MDAnalysis.analysis import distances
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

def get_contact_data(gro, xtc, system_name):
    """
    Computes residue-level contact occupancy (>10%) between TM2 (220-240) 
    and the surrounding environment (200-500).
    """
    if not os.path.exists(gro) or not os.path.exists(xtc):
        print(f"Error: Missing files for {system_name}")
        return pd.DataFrame()

    u = mda.Universe(gro, xtc)
    tm2 = u.select_atoms("resnum 220:240 and not name H*")
    env = u.select_atoms("resnum 200:500 and not name H*")
    
    tm2_resids = tm2.resnums
    env_resids = env.resnums
    unique_tm2 = np.unique(tm2_resids)
    unique_env = np.unique(env_resids)
    
    counts = np.zeros((len(unique_env), len(unique_tm2)))
    n_frames = len(u.trajectory)
    
    env_map = {rid: i for i, rid in enumerate(unique_env)}
    tm2_map = {rid: i for i, rid in enumerate(unique_tm2)}
    env_indices = np.array([env_map[r] for r in env_resids])
    tm2_indices = np.array([tm2_map[r] for r in tm2_resids])

    print(f"Processing {system_name} ({n_frames} frames)...")
    for ts in u.trajectory:
        dist = distances.distance_array(env.positions, tm2.positions)
        contacts = dist < 4.5
        
        env_c, tm2_c = np.where(contacts)
        if len(env_c) > 0:
            res_pairs = np.unique(np.column_stack((env_indices[env_c], tm2_indices[tm2_c])), axis=0)
            counts[res_pairs[:, 0], res_pairs[:, 1]] += 1

    occupancy = (counts / n_frames) * 100
    
    data = []
    for i, e_res in enumerate(unique_env):
        for j, t_res in enumerate(unique_tm2):
            if occupancy[i, j] > 10:
                data.append({
                    "System": system_name,
                    "TM2_Residue": t_res,
                    "Neighbor_Residue": e_res,
                    "Occupancy": occupancy[i, j]
                })
    return pd.DataFrame(data)

def plot_refined_heatmaps(df, output_path):
    """
    Generates balanced heatmaps using only interacting residues.
    Ensures identical grid sizes for both plots by using a dedicated colorbar axis.
    """
    if df.empty: return

    tm2_order = sorted(range(220, 241))
    neighbor_order = sorted(df['Neighbor_Residue'].unique())
    
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'DejaVu Sans'],
        'font.weight': 'bold',
        'axes.labelweight': 'bold',
        'axes.titleweight': 'bold'
    })

    # Reverted to previous stable size
    height_per_res = 1.0
    fig_height = max(16, len(neighbor_order) * height_per_res)
    
    from matplotlib.gridspec import GridSpec
    fig = plt.figure(figsize=(42, fig_height))
    gs = GridSpec(1, 3, width_ratios=[1, 1, 0.04], figure=fig)
    
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharey=ax1)
    cbar_ax = fig.add_subplot(gs[2])
    
    TITLE_SIZE = 60
    LABEL_SIZE = 52
    TICK_SIZE = 44

    systems = [('Apo', ax1), ('Bound', ax2)]
    
    for i, (sys_name, ax) in enumerate(systems):
        sys_df = df[df['System'] == sys_name]
        pivot = sys_df.pivot(index='Neighbor_Residue', columns='TM2_Residue', values='Occupancy')
        pivot = pivot.reindex(index=neighbor_order, columns=tm2_order)
        
        sns.heatmap(pivot, ax=ax, cmap="viridis", vmin=0, vmax=100,
                    cbar=(i == 1),
                    cbar_ax=cbar_ax if i == 1 else None,
                    cbar_kws={'label': 'Occupancy (%)'} if i == 1 else None,
                    linewidths=1.0, linecolor='lightgray', mask=pivot.isnull())
        
        ax.set_title(f"{sys_name} SPPL2a", fontsize=TITLE_SIZE, pad=60)
        ax.set_xlabel("TM2 Residue", fontsize=LABEL_SIZE, labelpad=30)
        
        if i == 0:
            ax.set_ylabel("Neighboring Residue", fontsize=LABEL_SIZE, labelpad=30)
            ax.set_yticklabels(neighbor_order, rotation=0, fontsize=TICK_SIZE)
        else:
            ax.set_ylabel("")
            plt.setp(ax.get_yticklabels(), visible=False)
            
        ax.set_xticklabels(tm2_order, rotation=45, fontsize=TICK_SIZE)

        for _, spine in ax.spines.items():
            spine.set_visible(True)
            spine.set_linewidth(5)

    cbar_ax.set_ylabel('Occupancy (%)', size=LABEL_SIZE, labelpad=40, weight='bold')
    cbar_ax.tick_params(labelsize=TICK_SIZE)

    plt.suptitle("TM2 Gating Contact Stability Fingerprint", fontsize=TITLE_SIZE + 15, y=0.992)
    plt.subplots_adjust(left=0.08, right=0.92, wspace=0.15, top=0.91)
    plt.savefig(output_path, bbox_inches='tight')
    tiff_path = output_path.replace('.pdf', '.tiff')
    plt.savefig(tiff_path, dpi=300, bbox_inches='tight')
    print(f'Heatmap also saved to {tiff_path}')
    print(f"Heatmap saved to {output_path}")

if __name__ == "__main__":
    os.makedirs('analysis/plots', exist_ok=True)
    df_apo = get_contact_data('9K92_main/step7_1.gro', 'analysis/9K92/reduced.xtc', 'Apo')
    df_bound = get_contact_data('9K93_main/step7_1.gro', 'analysis/9K93/reduced.xtc', 'Bound')
    if not df_apo.empty or not df_bound.empty:
        if not df_apo.empty:
            df_apo.to_csv('analysis/tm2_occupancy_apo.csv', index=False)
        if not df_bound.empty:
            df_bound.to_csv('analysis/tm2_occupancy_bound.csv', index=False)
        combined = pd.concat([df_apo, df_bound])
        combined.to_csv('analysis/tm2_occupancy.csv', index=False)
        plot_refined_heatmaps(combined, 'analysis/plots/tm2_contact_occupancy.pdf')
