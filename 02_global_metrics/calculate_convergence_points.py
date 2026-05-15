import numpy as np
import os

def read_xvg(filename):
    time, data = [], []
    if not os.path.exists(filename): return None, None
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith(('#', '@')): continue
            cols = line.strip().split()
            if len(cols) >= 2:
                time.append(float(cols[0]))
                data.append(float(cols[1]))
    return np.array(time), np.array(data)

def get_slope(x, y):
    n = len(x)
    if n < 2: return 0
    return (n * np.sum(x*y) - np.sum(x)*np.sum(y)) / (n * np.sum(x**2) - (np.sum(x)**2))

def find_strict_convergence(time, data, window_ns=100, slope_threshold=0.0001):
    dt = time[1] - time[0]
    window_steps = int(window_ns / dt)
    n_steps = len(data)
    
    # We want to find the first time T such that for EVERY sub-window from T to 1000, 
    # the slope is low. This prevents "accidental" early zero-slopes.
    
    for i in range(0, n_steps - window_steps):
        is_stable_from_here = True
        
        # Check all possible windows starting from index 'i' to the end
        for j in range(i, n_steps - window_steps, 50): # check every 50 steps for speed
            w_t = time[j : j + window_steps]
            w_d = data[j : j + window_steps]
            
            if abs(get_slope(w_t, w_d)) > slope_threshold:
                is_stable_from_here = False
                break
        
        if is_stable_from_here:
            # Final check on the total remaining segment
            if abs(get_slope(time[i:], data[i:])) < slope_threshold:
                return time[i]
                
    return None

systems = {
    "Apo SPPL2a (9K92)": ("analysis/9K92", "rmsd.xvg"),
    "Bound SPPL2a (9K93)": ("analysis/9K93", "rmsd.xvg"),
    "9K95 Complex": ("analysis/9K95", "rmsd.xvg"),
    "PS1 Region (9K95)": ("analysis/9K95", "rmsd_ps1_fixed.xvg")
}

print(f"{'System':<20} | {'Convergence Time (ns)':<22} | {'Final Mean RMSD (nm)'}")
print("-" * 75)

for name, (path, fname) in systems.items():
    t, d = read_xvg(os.path.join(path, fname))
    if d is None: continue
    
    # Threshold 0.0001 nm/ns is very strict (0.01 nm drift per 100 ns)
    conv_t = find_strict_convergence(t, d, window_ns=100, slope_threshold=0.0001)
    
    if conv_t is not None:
        mask = t >= conv_t
        final_mean = np.mean(d[mask])
        print(f"{name:<20} | {conv_t:>10.1f} ns            | {final_mean:.3f}")
    else:
        # Relax slightly if no point is found
        conv_t = find_strict_convergence(t, d, window_ns=100, slope_threshold=0.0002)
        if conv_t is not None:
             mask = t >= conv_t
             final_mean = np.mean(d[mask])
             print(f"{name:<20} | {conv_t:>10.1f} ns (relaxed)  | {final_mean:.3f}")
        else:
             print(f"{name:<20} | Did not reach plateau      | N/A")
