import pandas as pd
import numpy as np

def tidy_dssp(dssp_file, output_csv, offset, start_res=350, end_res=400):
    with open(dssp_file, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]
    
    data = []
    # Only take every 5th frame for the heatmap to keep file size manageable
    for i in range(0, len(lines), 5):
        line = lines[i]
        for j, char in enumerate(line):
            res_num = j + offset
            if start_res <= res_num <= end_res:
                data.append({
                    "Time": i, # Frame index as proxy for time
                    "Residue": res_num,
                    "Structure": char
                })
    
    df = pd.DataFrame(data)
    df.to_csv(output_csv, index=False)
    print(f"Tidy DSSP saved to {output_csv}")

# Corrected offsets: 9K92=169, 9K93=168
tidy_dssp("analysis/9K92/dssp.dat", "analysis/9K92/dssp_tidy.csv", offset=169)
tidy_dssp("analysis/9K93/dssp.dat", "analysis/9K93/dssp_tidy.csv", offset=168)