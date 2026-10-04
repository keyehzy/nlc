// Canonical labelling of small vertex-coloured graphs (at most 32 vertices).
//
// Two clusters with isomorphic structure have identical spectra and identical
// subcluster decompositions, so every cluster quantity is computed once per
// isomorphism class.  Block clusters are labelled through their spin-block
// incidence graph, with spins and blocks in different colours.
#pragma once

#include <array>
#include <compare>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace nlc {

inline constexpr int kMaxVertices = 32;
inline constexpr int kCertWords = (kMaxVertices * (kMaxVertices - 1) / 2 + 63) / 64;

using Edge = std::pair<int, int>;

// Initial colours in canonical order, then the upper-triangular adjacency bit
// string under the canonical labelling, pair (0,1) most significant.
struct Cert {
  int n = 0;
  std::array<std::uint8_t, kMaxVertices> colors{};
  std::array<std::uint64_t, kCertWords> bits{};
  auto operator<=>(const Cert&) const = default;
};

struct CertHash {
  std::size_t operator()(const Cert& c) const;
};

struct Canon {
  Cert cert;
  std::vector<int> perm;  // perm[v] = canonical index of vertex v
};

// Individualisation-refinement over the full search tree; the certificate is
// the minimum over all leaves, so it is a complete isomorphism invariant.
// `colors` (default: all zero) must be small non-negative integers; canonical
// indices list the vertices in increasing colour.
Canon canonical_form(int n, const std::vector<Edge>& edges, const std::vector<int>& colors = {});

// Disjoint pairs of twins (vertices u, v with the same neighbours apart from
// each other).  Swapping a twin pair is a graph automorphism, and swaps of
// disjoint pairs commute, so the pairs generate an abelian symmetry group.
std::vector<Edge> twin_pairs(int n, const std::vector<Edge>& edges);

}  // namespace nlc
