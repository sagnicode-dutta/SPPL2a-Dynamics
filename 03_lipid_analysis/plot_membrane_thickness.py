import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import os

# Enable Illustrator Editability (Type 42 fonts)
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# Global Style Settings (Matching the R theme_classic)
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Liberation Sans', 'DejaVu Sans']
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False

def plot_thickness():
    path = "analysis/membrane_thickness_detailed.csv"
    if not os.path.exists(path):
        print("Data file not found!")
        return

    df = pd.read_csv(path)
    
    # Use raw ANGSTROM values (do not divide by 10)
    # NORMALIZE: Divide Total by 2 to match Leaflet scale (~20 A)
    df["Total"] = df["Total"] / 2.0
        
    # Reshape for seaborn (melt) to match R's pivot_longer
    df_melt = df.melt(id_vars=["System"], value_vars=["Upper", "Lower", "Total"], 
                      var_name="Leaflet", value_name="Thickness_A")
    
    # Map system names and categories
    label_map = {
        "9K92": "Apo SPPL2a",
        "9K93": "Bound SPPL2a",
        "9K95": "Gamma-Secretase"
    }
    leaflet_map = {
        "Upper": "Upper leaflet",
        "Lower": "Lower leaflet",
        "Total": "Membrane"
    }
    
    df_melt['System_Label'] = df_melt['System'].map(label_map)
    df_melt['Leaflet_Label'] = df_melt['Leaflet'].map(leaflet_map)
    
    # Factor Order
    sys_order = ["Apo SPPL2a", "Bound SPPL2a", "Gamma-Secretase"]
    leaf_order = ["Upper leaflet", "Lower leaflet", "Membrane"]
    
    # Colors from reference script
    leaflet_colors = {
        "Membrane": "#6D2529",
        "Lower leaflet": "#2E5A87",
        "Upper leaflet": "#EB7822"
    }
    
    # A4 Dimensions: 11.69 x 8.27 (Landscape)
    fig, ax = plt.subplots(figsize=(11.69, 8.27))
    
    # Create Grouped Boxplot: X = System, Y = Thickness (in A)
    sns.boxplot(x='System_Label', y='Thickness_A', hue='Leaflet_Label', data=df_melt, 
                hue_order=leaf_order, order=sys_order, palette=leaflet_colors,
                width=0.7, linewidth=2, fliersize=4, 
                flierprops={"marker": 'o', "alpha": 0.5}, ax=ax)

    # Styling
    ax.set_title("Membrane Thickness", fontsize=32, pad=60)
    ax.set_xlabel("", fontsize=24)
    ax.set_ylabel("Thickness / Half-Thickness (Å)", fontsize=22)
    
    # Tighten Y-axis to focus on the data range (approx 18.5 to 20.5 A)
    ax.set_ylim(18.5, 20.5)
    
    ax.tick_params(labelsize=18, width=1.5, length=8)
    
    # Legend on top, pushed lower
    ax.legend(title="", loc='lower center', bbox_to_anchor=(0.5, 1.0), 
              ncol=3, fontsize=18, frameon=False)
    
    # Dashed Y-grid
    ax.grid(True, axis='y', color='grey', linestyle='--', alpha=0.3, zorder=0)
    ax.set_axisbelow(True)

    plt.tight_layout()
    plt.savefig("analysis/plots/membrane_thickness_boxplot.pdf", dpi=600, bbox_inches='tight')
    print("Membrane thickness boxplot (Å) saved to analysis/plots/membrane_thickness_boxplot.pdf")

if __name__ == "__main__":
    if not os.path.exists("analysis/plots"):
        os.makedirs("analysis/plots")
    plot_thickness()
