#!/bin/bash

# Configuration
SYSTEMS=("9k92" "9k93" "9k95")
# Times in ps (200ns, 400ns, 600ns, 800ns, 1000ns)
TIMES=(200000 400000 600000 800000 1000000)

for sys in "${SYSTEMS[@]}"; do
    echo "Processing $sys..."
    
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
    
    # 1. Clustering (Representative Structure)
    # Selection: TM_Core for fit, TM_Core for RMSD, Protein for output
    echo "Running clustering for $sys using reduced trajectory..."
    # Using specific group names and the -g flag for output to be explicit
    printf "TM_Core\nTM_Core\nProtein\n" | gmx cluster -s $tpr -f $xtc -n $ndx -cutoff 0.2 -method gromos \
        -cl analysis/snapshots/${sys}_all_clusters.pdb \
        -dist analysis/snapshots/${sys}_rmsd_dist.xvg \
        -sz analysis/snapshots/${sys}_cluster_size.xvg
    
    # Create a temporary protein-only structure for mapping
    printf "Protein\n" | gmx editconf -f $tpr -o analysis/snapshots/tmp_prot.gro
    
    # Extract only the first cluster (the centroid of the largest cluster)
    # The clusters file contains Protein atoms, and tmp_prot.gro also contains Protein atoms.
    # We still need to select the group for output in trjconv.
    echo "Protein" | gmx trjconv -f analysis/snapshots/${sys}_all_clusters.pdb -s analysis/snapshots/tmp_prot.gro -o analysis/snapshots/${sys}_representative.pdb -dump 0
    rm -f analysis/snapshots/tmp_prot.gro
    
    # 2. Extract specific snapshots from the reduced trajectory
    for t in "${TIMES[@]}"; do
        ns=$((t / 1000))
        echo "Extracting ${ns}ns snapshot for $sys..."
        # Selection: Protein for centering, Protein for output
        echo "Protein Protein" | gmx trjconv -s $tpr -f $xtc -dump $t -o analysis/snapshots/${sys}_${ns}ns.pdb -pbc mol -center -ur compact
    done
done

# Cleanup all-clusters file
rm -f analysis/snapshots/*_all_clusters.pdb
