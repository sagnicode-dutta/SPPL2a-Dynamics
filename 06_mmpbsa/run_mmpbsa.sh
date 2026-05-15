#!/bin/bash
# =============================================================================
# run_mmpbsa.sh - Template for Binding Energy Calculation (gmx_MMPBSA)
# =============================================================================
# This script generates the input file and runs the MM/PBSA calculation
# using the gmx_MMPBSA tool.
# =============================================================================

# ── 1. Write the Input Configuration File ──────────────────────────────
cat > mmpbsa.in << 'EOF'
&general
  sys_name    = "PROJECT_NAME",
  startframe  = 1,
  endframe    = 200,
  interval    = 1,
  verbose     = 2,
  keep_files  = 2,
/
&pb
  istrng      = 0.150,
  fillratio   = 4.0,
  indi        = 4.0,
  inp         = 2,
/
&decomp
  idecomp     = 1,
  dec_verbose = 3,
  print_res   = "within 6",
/
EOF

# ── 2. Run the Calculation ────────────────────────────────────────────
# Note: Update file names below to match your system files.
# -cs: complex.tpr
# -ct: trajectory.xtc
# -cp: topology.top
# -ci: index.ndx
# -cg: receptor_group ligand_group

gmx_MMPBSA -O \
  -i mmpbsa.in \
  -cs complex.tpr \
  -ct trajectory.xtc \
  -cp complex.top \
  -ci index.ndx \
  -cg 1 13 \
  -o FINAL_RESULTS.dat \
  -do FINAL_DECOMP.dat \
  -nogui
