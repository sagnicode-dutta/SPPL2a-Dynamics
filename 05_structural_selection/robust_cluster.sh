#!/bin/bash
set -euo pipefail

# Configuration
SYSTEMS=("9k92" "9k93" "9k95")
CUTOFF=0.20
METHOD="gromos"
FIT_GROUP="TM_Core"
OUTPUT_GROUP="Protein"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
error() { echo "[ERROR] $*" >&2; exit 1; }

for sys in "${SYSTEMS[@]}"; do
    log "Processing $sys..."
    
    # Path setup
    if [ "$sys" == "9k92" ]; then
        xtc="analysis/9K92/reduced.xtc"
        tpr="9K92_main/step7_1.tpr"
    elif [ "$sys" == "9k93" ]; then
        xtc="analysis/9K93/reduced.xtc"
        tpr="9K93_main/step7_1.tpr"
    elif [ "$sys" == "9k95" ]; then
        xtc="analysis/9K95/reduced.xtc"
        tpr="9K95_main/step7_1.tpr"
    fi
    
    ndx="analysis/snapshots/${sys}_core.ndx"
    out_pdb="analysis/snapshots/${sys}_all_clusters.pdb"
    rep_pdb="analysis/snapshots/${sys}_representative.pdb"
    
    log "Step 1: Running gmx cluster for $sys (RMSD Matrix)..."
    # Note: GROMACS 2025.2 combined Fit/RMSD prompt. Providing TM_Core then Protein.
    # If it asks for 3 groups, this might need adjustment, but usually it's 2 in this version.
    echo "TM_Core Protein" | gmx cluster -s "$tpr" -f "$xtc" -n "$ndx" -cutoff "$CUTOFF" -method "$METHOD" -cl "$out_pdb" -g "analysis/snapshots/${sys}_cluster.log" -dist "analysis/snapshots/${sys}_dist.xvg" -sz "analysis/snapshots/${sys}_size.xvg" -quiet

    log "Step 2: Extracting centroid using awk..."
    if [[ -f "$out_pdb" ]]; then
        awk '/^MODEL        1/,/^ENDMDL/' "$out_pdb" | grep -v "^MODEL\|^ENDMDL" > "$rep_pdb"
        log "✓ Success: $rep_pdb generated."
    else
        error "Cluster output $out_pdb was not created for $sys"
    fi
done

log "All systems processed successfully."
