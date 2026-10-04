# nlc — numerical linked-cluster expansion in C++

A C++20 reproduction of M. Rigol, T. Bryant and R. R. P. Singh,
*Numerical Linked-Cluster Approach to Quantum Lattice Models*,
PRL **97**, 187202 (2006): finite-temperature energy, entropy and specific
heat of spin-1/2 models in the thermodynamic limit from full exact
diagonalisation of finite clusters,

    H = sum_<ij> [ S^z_i S^z_j + J_perp (S^x_i S^x_j + S^y_i S^y_j) ] - h_x sum_i S^x_i - h_z sum_i S^z_i

on the kagome, triangular and square lattices.

## Build and run

Requires CMake >= 3.24, a C++20 compiler, Boost headers and LAPACK
(Accelerate on macOS; `brew install cmake ninja boost`).

    cmake -S . -B build -G Ninja
    cmake --build build
    ctest --test-dir build          # unit tests (a few seconds)
    tools/run_all.sh                # data/ for Figs. 1-5, then figures/

Individual runs:

    ./build/nlc clusters --lattice kagome --order 8
    ./build/nlc run --lattice kagome --order 8 --jperp 0,0.25,1 --out data/xxz
    ./build/nlc run --lattice kagome --order 6 --hx 0.5 --hz 0.25 --out data/tlfim
    ./build/nlc run --lattice square --order 13 --jperp 1 --out data/sq

`--jperp`, `--hx`, `--hz` take comma-separated lists and every combination
is written to `PREFIX_jp<J>_hx<hx>_hz<hz>.csv`.  Columns: `T`, then for
`E`, `S`, `C` (per site) the bare partial sums through each order
(`S_6` = all clusters with up to 6 triangles), Wynn's epsilon algorithm with
c cycles (`S_wynn2`), Brezinski's theta algorithm (`S_brez2`) and the Euler
transformation from order `--euler-start` (default 2) on (`S_euler`).
Other options: `--tmin --tmax --nt` (log grid, default 1e-3..100, 401
points), `--threads`, `--memory-gb` (cap on concurrently held dense blocks,
default 8).

## Method

* **Clusters** (`lattice`, `expansion`).  Kagome and triangular use the
  triangle expansion of the paper: clusters are connected sets of triangles
  (all triangles of the kagome lattice, the up-pointing ones of the
  triangular lattice; either way each bond lies in exactly one triangle),
  their Hamiltonian holds only the bonds of those triangles, and the single
  site is order 0.  The triangles are enumerated as site clusters of the
  lattice of triangle centres (honeycomb for kagome, a triangular lattice
  with the C3v point group of an up-triangle for the triangular lattice),
  modulo translations and point group.  The square lattice uses the site
  expansion (connected site sets with all bonds between them).  Clusters
  are grouped into topological classes by an individualisation-refinement
  canonical form of the spin–triangle incidence graph (or of the bond graph
  for site clusters); every proper connected subcluster of a class is
  enumerated and classified once, giving the inclusion–exclusion
  coefficients of Eq. (2).
* **Exact diagonalisation** (`ed`).  Each class is block-diagonalised by
  total S^z (when h_x = 0), the global spin flip (between ±S^z, or as a parity
  in the S^z = 0 / h_x != 0 space, when h_z = 0), and the parities of
  disjoint twin swaps of the cluster graph (e.g. the two free corners of a
  dangling triangle), which form an abelian group.  Blocks are diagonalised
  densely with LAPACK's two-stage `dsyevd_2stage` (eigenvalues only): the
  one-stage reduction streams the whole matrix through memory for every
  column and stops scaling once several blocks run concurrently.
* **Thermodynamics and sums** (`thermo`, `nlc`, `resum`).  E, S and C of
  each cluster are computed from its spectrum, measured from the ground
  state so that T = 10^-3 is accurate; weights are subtracted in class
  order and summed with the lattice constants order by order.

Cluster counts (classes / embeddings per unit cell):

| order | kagome triangles | triangular up-triangles | square sites |
|------:|-----------------:|------------------------:|-------------:|
| 1 | 1 / 2 | 1 / 1 | 1 / 1 |
| 2 | 1 / 3 | 1 / 3 | 1 / 2 |
| 3 | 1 / 6 | 3 / 11 | 1 / 6 |
| 4 | 2 / 14 | 5 / 44 | 3 / 19 |
| 5 | 2 / 36 | 12 / 186 | 4 / 63 |
| 6 | 5 / 94 | 35 / 814 | 10 / 216 |
| 7 | 7 / 250 | 98 / 3652 | 19 / 760 |
| 8 | 15 / 675 | 299 / 16689 | 51 / 2725 |
| 13 | | | 5450 / 1903890 |

The kagome class counts are the paper's 1, 1, 1, 2, 2, 5, 7, 15.

## Results

`tools/run_all.sh` produces `figures/fig1..fig5` (Apple M4, 10 cores):

| figure | model | expansion | time |
|---|---|---|---|
| 1 | kagome Ising, h_x = 0 ... 1 | 6 triangles | 3 s |
| 2, 3 | kagome XXZ, J_perp = 0 ... 1 | 8 triangles (17 spins) | 8 min |
| 4 | kagome Ising, h_x = 0.5, h_z = 0.25 | 6 triangles | 8 s |
| 5 | triangular / square Heisenberg | 8 up-triangles / 13 sites | 9 / 5 min |

Checks against numbers quoted in the paper:

* kagome Ising antiferromagnet, S(T = 10^-3): 0.50136 through 5 triangles
  (Pauling) and 0.50182 with the 12-site hexagon (exact 0.50183); both are
  unit tests.
* Fig. 1: for h_x = 0, 0.5, 1 the 5- and 6-triangle sums coincide at all T;
  h_x = 0.25 shows the double-peaked C_v; h_x = 0.01 the two separated peaks
  with the low-T one beyond convergence below T ~ 0.005.
* Figs. 2-3: the kagome Heisenberg C_v peak is 0.189 at T = 0.68 and the 7-
  and 8-triangle sums agree down to T ~ 0.4, with the 8-triangle sum dipping
  to 0.13 at T ~ 0.3 before it diverges.
* Fig. 4: bare sums of 3-6 triangles diverge below T ~ 1 while Wynn_2 and
  Euler agree at all temperatures (peak C_v = 0.26 at T ~ 0.25).
* Fig. 5: square-lattice Wynn_2 / Brezinski_2 agree down to T ~ 0.3
  (C_v peak 0.456 at T ~ 0.58); triangular Wynn_3 and Wynn_4 agree down to
  T ~ 0.25, below which Wynn_4 levels at S ~ 0.2 and Wynn_3 turns up.

Not reproduced: the HTE series and Pade approximants of Fig. 3 (from Elstner
and Young, PRB 50, 6871) and the Bernu-Misguich curves of Fig. 5, which are
external data, and the schematic phase-diagram inset of Fig. 4.

## Code from nlce_cpp

`graph` (canonical labelling) and `lattice` (lattices, point groups, site
cluster enumeration) are adapted from `nlce_cpp`; the labelling was extended
to vertex colours and 32 vertices.  Everything else is new.

## License

MIT; see `LICENSE`.
