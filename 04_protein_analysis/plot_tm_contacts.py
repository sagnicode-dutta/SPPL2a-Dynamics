import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

def plot_interface(df_apo, df_bound, interface, output_dir, title):
    # Pivot to matrix
    m_apo = df_apo[df_apo['Interface'] == interface].pivot(index='Residue_A', columns='Residue_B', values='Occupancy').fillna(0)
    m_bound = df_bound[df_bound['Interface'] == interface].pivot(index='Residue_A', columns='Residue_B', values='Occupancy').fillna(0)
    
    # Ensure both matrices have same residues
    all_a = sorted(list(set(m_apo.index) | set(m_bound.index)))
    all_b = sorted(list(set(m_apo.columns) | set(m_bound.columns)))
    
    m_apo = m_apo.reindex(index=all_a, columns=all_b).fillna(0)
    m_bound = m_bound.reindex(index=all_a, columns=all_b).fillna(0)
    
    delta = m_bound - m_apo
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    sns.heatmap(m_apo, ax=axes[0], cmap="YlGnBu", vmin=0, vmax=1)
    axes[0].set_title(f"{title} (Apo)")
    
    sns.heatmap(m_bound, ax=axes[1], cmap="YlGnBu", vmin=0, vmax=1)
    axes[1].set_title(f"{title} (Bound)")
    
    sns.heatmap(delta, ax=axes[2], cmap="RdBu_r", vmin=-0.5, vmax=0.5)
    axes[2].set_title("Delta (Bound - Apo)")
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/contacts_{interface.lower()}_{title.replace(' ', '_')}.pdf")
    plt.close()

# Load data
df_9k92 = pd.read_csv('analysis/9K92/tm_contacts.csv')
df_9k93 = pd.read_csv('analysis/9K93/tm_contacts.csv')
df_9k95 = pd.read_csv('analysis/9K95/tm_contacts.csv')

# Plot SPPL2a differences
for interface in ["TM3-TM4", "TM6-TM7", "TM8-TM9"]:
    plot_interface(df_9k92, df_9k93, interface, "analysis/plots", "SPPL2a")

# PS1 vs Bound SPPL2a is harder to plot as delta because residues differ.
# We will just plot them side by side.

def plot_side_by_side(df_sppl, df_ps1, interface, output_dir):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    m_sppl = df_sppl[df_sppl['Interface'] == interface].pivot(index='Residue_A', columns='Residue_B', values='Occupancy').fillna(0)
    m_ps1 = df_ps1[df_ps1['Interface'] == interface].pivot(index='Residue_A', columns='Residue_B', values='Occupancy').fillna(0)
    
    sns.heatmap(m_sppl, ax=axes[0], cmap="YlGnBu", vmin=0, vmax=1)
    axes[0].set_title(f"SPPL2a (Bound) {interface}")
    
    sns.heatmap(m_ps1, ax=axes[1], cmap="YlGnBu", vmin=0, vmax=1)
    axes[1].set_title(f"PS1 (Bound) {interface}")
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/contacts_comparison_{interface.lower()}.pdf")
    plt.close()

for interface in ["TM3-TM4", "TM6-TM7", "TM8-TM9"]:
    plot_side_by_side(df_9k93, df_9k95, interface, "analysis/plots")
