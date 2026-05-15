import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import os

# Enable Illustrator Editability (Type 42 fonts)
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# Global Style Settings for High Visibility
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Liberation Sans', 'DejaVu Sans']
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'
plt.rcParams['lines.linewidth'] = 4
plt.rcParams['axes.linewidth'] = 3
plt.rcParams['xtick.major.width'] = 3
plt.rcParams['ytick.major.width'] = 3
plt.rcParams['xtick.major.size'] = 10
plt.rcParams['ytick.major.size'] = 10

def plot_scd():
    systems = [
        {"name": "9K92", "label": "Apo SPPL2a", "color": "#FF0046"},
        {"name": "9K93", "label": "Bound SPPL2a", "color": "#53C719"},
        {"name": "9K95", "label": "Gamma-Secretase", "color": "#9929EA"},
    ]

    # Horizontal Layout (1x2) as requested
    # Using A4 Landscape dimensions (11.69 x 8.27) to maintain A4 area
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))

    # Plot sn-1 (Palmitoyl)
    for sys_info in systems:
        path = f"analysis/{sys_info['name']}_scd_palmitoyl.csv"
        if os.path.exists(path):
            df = pd.read_csv(path)
            ax1.plot(df['Carbon'], df['SCD'], 'o-', color=sys_info['color'], 
                     label=sys_info['label'], markersize=14, markeredgecolor='black', markeredgewidth=2)

    ax1.set_title("sn-1 (Palmitoyl Chain)", fontsize=32, pad=20)
    ax1.set_xlabel("Carbon Number", fontsize=28)
    ax1.set_ylabel("$|S_{CD}|$", fontsize=28)
    ax1.set_ylim(0, 0.25)
    ax1.tick_params(labelsize=22)
    ax1.grid(True, linestyle='--', alpha=0.3)
    ax1.legend(fontsize=20)

    # Plot sn-2 (Oleoyl)
    for sys_info in systems:
        path = f"analysis/{sys_info['name']}_scd_oleoyl.csv"
        if os.path.exists(path):
            df = pd.read_csv(path)
            ax2.plot(df['Carbon'], df['SCD'], 's-', color=sys_info['color'], 
                     label=sys_info['label'], markersize=14, markeredgecolor='black', markeredgewidth=2)

    ax2.set_title("sn-2 (Oleoyl Chain)", fontsize=32, pad=20)
    ax2.set_xlabel("Carbon Number", fontsize=28)
    ax2.set_ylabel("$|S_{CD}|$", fontsize=28)
    ax2.set_ylim(0, 0.25)
    ax2.tick_params(labelsize=22)
    ax2.grid(True, linestyle='--', alpha=0.3)
    ax2.legend(fontsize=20)

    plt.tight_layout(pad=4.0)
    plt.savefig("analysis/plots/lipid_order_parameter.pdf", dpi=600, bbox_inches='tight')
    print("Lipid order parameter plot saved to analysis/plots/lipid_order_parameter.pdf")

if __name__ == "__main__":
    if not os.path.exists("analysis/plots"):
        os.makedirs("analysis/plots")
    plot_scd()
