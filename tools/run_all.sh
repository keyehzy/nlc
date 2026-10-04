#!/bin/sh
# Generates the data for Figs. 1-5 into data/ and plots them into figures/.
# Run from the project root after building (cmake --build build).
set -e
NLC=./build/nlc
mkdir -p data
GRID="--tmin 1e-3 --tmax 100 --nt 501"

# Fig. 1: transverse-field Ising model on the kagome lattice, 5 and 6 triangles.
$NLC run --lattice kagome --order 6 --hx 0,0.01,0.25,0.5,1.0 $GRID --out data/fig1_kagome_tfim
# Figs. 2 and 3: XXZ models on the kagome lattice, 7 and 8 triangles.
$NLC run --lattice kagome --order 8 --jperp 0,0.01,0.1,0.25,1 $GRID --out data/fig2_kagome_xxz
# Fig. 4: Ising model in hx = 0.5 and hz = 0.25 on the kagome lattice.
$NLC run --lattice kagome --order 6 --hx 0.5 --hz 0.25 $GRID --out data/fig4_kagome_tlfim
# Fig. 5: Heisenberg models, 8 up-triangles (triangular) and 13 sites (square).
$NLC run --lattice triangular --order 8 --jperp 1 $GRID --out data/fig5_triangular
$NLC run --lattice square --order 13 --jperp 1 $GRID --out data/fig5_square

python3 tools/plot_figures.py --data data --out figures
