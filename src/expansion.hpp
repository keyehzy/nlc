// Linked-cluster expansions: topological cluster classes, their lattice
// constants and their subcluster decompositions.
//
// P(L)/N = sum_c L(c) W(c),  W(c) = P(c) - sum_{s subset c} W(s)   (Eqs. 1, 2)
//
// Two families of clusters are supported.
//
// * Triangle expansions (kagome, triangular).  A cluster is a connected set
//   of triangles (every triangle of the kagome lattice, the up-pointing ones
//   of the triangular lattice; either way each bond lies in exactly one
//   triangle) and its Hamiltonian holds the bonds of those triangles only.
//   The single site is a cluster of order 0 and the order of any other
//   cluster is its number of triangles.  Clusters are classified by their
//   spin-triangle incidence graph, which fixes both the Hamiltonian and the
//   set of subclusters.
// * Site expansion (square).  A cluster is a connected set of sites with
//   every lattice bond between them, of order equal to its number of sites.
#pragma once

#include "graph.hpp"

#include <string>
#include <utility>
#include <vector>

namespace nlc {

struct ClusterClass {
  int order = 0;
  int nspins = 0;
  std::vector<Edge> bonds;                  // spin graph, canonical labels
  std::vector<std::vector<int>> blocks;     // triangles (triangle expansions only)
  double lattice_constant = 0.0;            // embeddings per lattice site
  std::vector<std::pair<int, int>> subclusters;  // (class, count) over proper subclusters
};

struct Expansion {
  std::string lattice;
  int max_order = 0;
  std::vector<ClusterClass> classes;  // nondecreasing order; subclusters precede their parents
  std::vector<int> classes_per_order;
  std::vector<long long> embeddings_per_cell;  // [order], lattice clusters per unit cell
};

// lattice: "kagome", "triangular" (triangle expansions) or "square" (site
// expansion).  `threads` parallelises the subcluster enumeration.
Expansion build_expansion(const std::string& lattice, int max_order, int threads = 1);

// Same expansion enumerated with the trivial point group (cross-check only).
Expansion build_expansion_translations_only(const std::string& lattice, int max_order);

}  // namespace nlc
