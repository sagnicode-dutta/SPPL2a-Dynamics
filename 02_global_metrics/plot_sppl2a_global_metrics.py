import numpy as np
import matplotlib.pyplot as plt
import os

# Ultra-High-Visibility Global Settings (matching previous plots)
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Liberation Sans', 'DejaVu Sans']
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'
plt.rcParams['lines.linewidth'] = 4
plt.rcParams['axes.linewidth'] = 3

def read_xvg(filename):
    time, data = [], []
    if not os.path.exists(filename): return np.array([]), np.array([])
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith(('#', '@')): continue
            cols = line.strip().split()
            if len(cols) >= 2:
                time.append(float(cols[0]))
                data.append(float(cols[1]))
    return np.array(time), np.array(data)

# Colors
color_apo = '#FF0046'
color_bound = '#53C719'

# Data Paths
metrics = ['rmsd', 'gyrate', 'sasa', 'hbond']
titles = ['RMSD (Backbone)', 'Radius of Gyration', 'SASA (Total)', 'Internal H-bonds']
ylabels = ['RMSD (nm)', 'Rg (nm)', 'Area (nm$^2$)', 'Count']

fig, axes = plt.subplots(2, 2, figsize=(24, 18))
axes = axes.flatten()

for i, metric in enumerate(metrics):
    t_92, d_92 = read_xvg(f'analysis/9K92/{metric}.xvg')
    t_93, d_93 = read_xvg(f'analysis/9K93/{metric}.xvg')
    
    if len(d_92) > 0:
        axes[i].plot(t_92, d_92, color=color_apo, label='Apo SPPL2a (9K92)', alpha=0.8, lw=3)
    if len(d_93) > 0:
        axes[i].plot(t_93, d_93, color=color_bound, label='Bound SPPL2a (9K93)', alpha=0.8, lw=3)
    
    axes[i].set_title(titles[i], fontsize=32, pad=20)
    axes[i].set_xlabel('Time (ns)', fontsize=28)
    axes[i].set_ylabel(ylabels[i], fontsize=28)
    axes[i].tick_params(labelsize=24, length=10, width=3)
    
    if i == 0:
        axes[i].legend(fontsize=20, frameon=True, edgecolor='black')

plt.suptitle('Global Stability & Metrics: SPPL2a Apo vs Bound', fontsize=40, fontweight='bold', y=0.98)
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig('analysis/plots/sppl2a_global_metrics_comparison.pdf', dpi=300)
print("Saved analysis/plots/sppl2a_global_metrics_comparison.pdf")
