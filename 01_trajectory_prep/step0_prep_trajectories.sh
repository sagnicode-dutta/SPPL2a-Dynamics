#!/bin/bash
# =============================================================================
# step0_prep_trajectories.sh
# =============================================================================
# This script prepares the protein-only trajectory for analysis by removing 
# PBC effects, centering the protein, and fitting to the initial frame.
# =============================================================================
set -e

# Update BASE to your project directory
BASE="./simulations"
SYSTEMS=("9K92_main" "9K93_main" "9K95_main")

PROTEIN_GROUP=1      # Protein
CA_GROUP=3           # C-alpha

for SYS in "${SYSTEMS[@]}"; do
    echo "Preparing: $SYS"

    SYS_DIR="$BASE/$SYS"
    OUT="./analysis/$SYS"
    mkdir -p "$OUT"

    TPR="$SYS_DIR/step7_1.tpr"
    XTC="$SYS_DIR/step7_1.xtc"

    # 1. Protein-only trajectory (Remove PBC + Center)
    echo "  [1/2] Removing PBC and centering protein..."
    printf "${PROTEIN_GROUP}\n${PROTEIN_GROUP}\n" | gmx trjconv -s "$TPR" -f "$XTC" -o "$OUT/traj_center.xtc" -pbc mol -center &> /dev/null

    # 2. Fit to CA atoms and output full protein
    echo "  [2/2] Fitting to C-alpha and extracting protein-only trajectory..."
    printf "${CA_GROUP}\n${PROTEIN_GROUP}\n" | gmx trjconv -s "$TPR" -f "$OUT/traj_center.xtc" -o "$OUT/protein_only.xtc" -fit rot+trans &> /dev/null
    
    # 3. Reference Structure
    printf "${PROTEIN_GROUP}\n" | gmx trjconv -s "$TPR" -f "$OUT/protein_only.xtc" -o "$OUT/ref_protein.pdb" -dump 0 &> /dev/null

    rm -f "$OUT/traj_center.xtc"
    echo "  DONE: $SYS"
done
