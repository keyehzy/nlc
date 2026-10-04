// Full exact diagonalisation of spin-1/2 XXZ clusters in transverse and
// longitudinal fields,
//
//   H = sum_<ij> [ Jz S^z_i S^z_j + Jperp (S^x_i S^x_j + S^y_i S^y_j) ]
//       - hx sum_i S^x_i - hz sum_i S^z_i                              (Eq. 3)
//
// The Hilbert space is block-diagonalised by every symmetry that is cheap to
// exploit: total S^z when hx = 0, the global spin flip when hz = 0 (between
// S^z = +-m sectors, or as a parity at m = 0), and the parities of disjoint
// twin swaps of the cluster graph (the two free corners of a dangling
// triangle, for instance).  Each block is diagonalised densely with LAPACK.
#pragma once

#include "graph.hpp"

#include <cstddef>
#include <vector>

namespace nlc {

struct Model {
  double jz = 1.0;
  double jperp = 0.0;
  double hx = 0.0;
  double hz = 0.0;
};

// Distinct levels in increasing order with their degeneracies.
struct Spectrum {
  std::vector<double> energy;
  std::vector<double> degeneracy;
  std::size_t largest_block = 0;
};

Spectrum spectrum(int n, const std::vector<Edge>& bonds, const Model& model);

// The full 2^n Hamiltonian diagonalised without symmetries (tests only).
Spectrum spectrum_without_symmetries(int n, const std::vector<Edge>& bonds, const Model& model);

// Caps the memory held by dense blocks across threads; a block larger than
// the cap runs alone.  The default is 8 GiB.
void set_memory_budget(std::size_t bytes);

// Makes LAPACK single-threaded where that is configurable, so clusters can
// be diagonalised in parallel.
void use_single_threaded_lapack();

}  // namespace nlc
