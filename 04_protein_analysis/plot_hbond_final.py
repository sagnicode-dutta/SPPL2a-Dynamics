import numpy as np
import matplotlib.pyplot as plt
import os

# Ultra-High-Visibility Global Settings
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Liberation Sans', 'DejaVu Sans']
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'
plt.rcParams['lines.linewidth'] = 4
plt.rcParams['axes.linewidth'] = 3

def read_xvg(filename):
    time, data = [], []
    if not os.path.exists(filename): 
        print(f"Warning: {filename} not found.")
        return np.array([]), np.array([])
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith(('#', '@')): continue
            cols = line.strip().split()
            if len(cols) >= 2:
                time.append(float(cols[0]))
                data.append(float(cols[1]))
    return np.array(time), np.array(data)

def smooth_data(data, window=51):
    if len(data) < window: return data
    pad_width = window // 2
    padded_data = np.pad(data, pad_width, mode='reflect')
    smoothed = np.convolve(padded_data, np.ones(window)/window, mode='valid')
    return smoothed[:len(data)]

# Colors
color_apo = '#FF0046'     # Pink
color_bound = '#53C719'   # Green
color_ps1 = '#9929EA'     # Purple
color_complex = '#3C0266' # Dark Purple
WINDOW_SIZE = 51

# 1. Plot SPPL2a H-bonds (Apo vs Bound)
fig, ax = plt.subplots(figsize=(15, 10))
t_92, d_92 = read_xvg('analysis/9K92/hbond.xvg')
t_93, d_93 = read_xvg('analysis/9K93/hbond.xvg')

if len(d_92) > 0:
    ax.plot(t_92, d_92, color=color_apo, alpha=0.3, lw=2.0, label='_nolegend_')
    ax.plot(t_92, smooth_data(d_92, WINDOW_SIZE), color=color_apo, alpha=1.0, lw=4, label='Apo SPPL2a (9K92)')
if len(d_93) > 0:
    ax.plot(t_93, d_93, color=color_bound, alpha=0.3, lw=2.0, label='_nolegend_')
    ax.plot(t_93, smooth_data(d_93, WINDOW_SIZE), color=color_bound, alpha=1.0, lw=4, label='Bound SPPL2a (9K93)')

ax.set_title('SPPL2a Internal Hydrogen Bonds', fontsize=36, pad=50)
ax.set_xlabel('Time (ns)', fontsize=32, labelpad=15)
ax.set_ylabel('H-bond Count', fontsize=32, labelpad=15)
ax.set_xlim(0, 1000)
ax.set_ylim(220, 280)

ax.tick_params(labelsize=28, length=12, width=3)
ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=2, fontsize=24, frameon=False)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig('analysis/plots/sppl2a_hbond_final.pdf', dpi=300, bbox_inches='tight')
print("Saved analysis/plots/sppl2a_hbond_final.pdf")
plt.close()

# 2. Plot PS1 vs Complex H-bonds (Broken Axis)
def plot_broken_hbond(complex_file, ps1_file, title, ylabel, filename):
    fig, (ax_top, ax_bot) = plt.subplots(2, 1, sharex=True, figsize=(15, 10), gridspec_kw={'height_ratios': [1, 1]})
    
    # Create a big invisible axis for a centered Y-label
    big_ax = fig.add_subplot(111, frameon=False)
    big_ax.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    big_ax.set_ylabel(ylabel, fontsize=32, fontweight='bold', labelpad=80)
    
    plt.subplots_adjust(hspace=0.2, top=0.82)
    
    t_comp, d_comp = read_xvg(complex_file)
    t_ps1, d_ps1 = read_xvg(ps1_file)
    
    line1, = ax_top.plot([], [], color=color_complex, lw=4, label='Entire Complex (9K95)')
    line2, = ax_top.plot([], [], color=color_ps1, lw=4, label='PS1 Region')

    for ax in [ax_top, ax_bot]:
        if len(d_comp) > 0:
            ax.plot(t_comp, d_comp, color=color_complex, alpha=0.3, lw=2.0)
            ax.plot(t_comp, smooth_data(d_comp, WINDOW_SIZE), color=color_complex, alpha=1.0, lw=4)
        if len(d_ps1) > 0:
            ax.plot(t_ps1, d_ps1, color=color_ps1, alpha=0.3, lw=2.0)
            ax.plot(t_ps1, smooth_data(d_ps1, WINDOW_SIZE), color=color_ps1, alpha=1.0, lw=4)

    ax_bot.set_ylim(230, 290)
    ax_top.set_ylim(980, 1100)
    ax_bot.set_yticks([240, 260, 280])
    ax_top.set_yticks([1000, 1040, 1080])
    
    ax_top.spines['bottom'].set_visible(False); ax_top.spines['top'].set_visible(False); ax_top.spines['right'].set_visible(False)
    ax_bot.spines['top'].set_visible(False); ax_bot.spines['right'].set_visible(False)
    
    ax_top.tick_params(labelbottom=False, bottom=False, labelsize=28, length=12, width=3, pad=10)
    ax_bot.tick_params(labelsize=28, length=12, width=3, pad=10)
    ax_bot.set_xlim(0, 1000)
    ax_bot.set_xlabel('Time (ns)', fontsize=32, labelpad=15)
    
    fig.suptitle(title, fontsize=36, fontweight='bold', y=0.96)
    fig.legend(handles=[line1, line2], loc='upper center', bbox_to_anchor=(0.5, 0.91), ncol=2, fontsize=24, frameon=False)

    d = .015
    kwargs = dict(transform=ax_top.transAxes, color='black', clip_on=False, lw=3)
    ax_top.plot((-d, +d), (-d, +d), **kwargs)
    kwargs.update(transform=ax_bot.transAxes)
    ax_bot.plot((-d, +d), (1 - d, 1 + d), **kwargs)
    
    plt.savefig(f'analysis/plots/{filename}', dpi=300, bbox_inches='tight')
    plt.close()

plot_broken_hbond('analysis/9K95/hbond.xvg', 'analysis/9K95/hbond_ps1_fixed.xvg', 
                  'Gamma-Secretase H-bond Dynamics', 'H-bond Count', '9k95_hbond_final.pdf')
print("Saved analysis/plots/9k95_hbond_final.pdf")
