#!/bin/bash

# Define systems and their PDB IDs
declare -A system_names=( ["9K92"]="Apo SPPL2a" ["9K93"]="Bound SPPL2a" ["9K95"]="Bound PS1" )

for sys in 9K92 9K93 9K95; do
    echo "Running analysis for $sys (${system_names[$sys]})..."
    DIR="analysis/$sys"
    TPR="${sys}_main/step7_1.tpr"
    XTC="$DIR/reduced.xtc"

    # 1. RMSD (Backbone-Backbone) - Group 4 4
    echo "4 4" | gmx rms -s "$TPR" -f "$XTC" -o "$DIR/rmsd.xvg" -tu ns

    # 2. RMSF (C-alpha) - Group 3
    echo "3" | gmx rmsf -s "$TPR" -f "$XTC" -o "$DIR/rmsf.xvg" -res

    # 3. SASA (Protein) - Group 1
    echo "1" | gmx sasa -s "$TPR" -f "$XTC" -o "$DIR/sasa.xvg" -tu ns

    # 4. Radius of Gyration (Protein) - Group 1
    echo "1" | gmx gyrate -s "$TPR" -f "$XTC" -o "$DIR/gyrate.xvg" -tu ns

    # 5. RMS Dist (Backbone) - Group 4
    echo "4" | gmx rmsdist -s "$TPR" -f "$XTC" -o "$DIR/rmsdist.xvg"
done
