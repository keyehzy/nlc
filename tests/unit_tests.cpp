// Unit tests: labelling, cluster enumeration, exact diagonalisation,
// thermodynamics, extrapolation, and the T -> 0 entropy of the kagome Ising
// antiferromagnet quoted in the paper.
#include "ed.hpp"
#include "expansion.hpp"
#include "graph.hpp"
#include "nlc.hpp"
#include "resum.hpp"
#include "thermo.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <numeric>
#include <random>
#include <string>
#include <vector>

using namespace nlc;

namespace {

int failures = 0;

void check(bool ok, const std::string& what) {
  if (!ok) {
    ++failures;
    std::fprintf(stderr, "FAIL: %s\n", what.c_str());
  }
}

bool close(double a, double b, double tol) { return std::abs(a - b) <= tol; }

void test_canonical_form() {
  std::mt19937 rng(7);
  for (int trial = 0; trial < 200; ++trial) {
    const int n = 3 + trial % 14;
    std::vector<Edge> edges;
    std::vector<int> colors(n);
    for (int i = 0; i < n; ++i) {
      colors[i] = rng() % 2;
      for (int j = i + 1; j < n; ++j) {
        if (rng() % 3 == 0) edges.emplace_back(i, j);
      }
    }
    std::vector<int> perm(n);
    std::iota(perm.begin(), perm.end(), 0);
    std::shuffle(perm.begin(), perm.end(), rng);
    std::vector<Edge> moved;
    for (auto [u, v] : edges) moved.emplace_back(perm[u], perm[v]);
    std::vector<int> moved_colors(n);
    for (int i = 0; i < n; ++i) moved_colors[perm[i]] = colors[i];
    check(canonical_form(n, edges, colors).cert == canonical_form(n, moved, moved_colors).cert,
          "certificate invariant under relabelling");
    if (!edges.empty()) {
      auto fewer = edges;
      fewer.pop_back();
      check(!(canonical_form(n, edges).cert == canonical_form(n, fewer).cert), "certificate separates graphs");
    }
  }
  // A path a-b-c with colours (1,0,0) and (0,1,0) are not isomorphic.
  const std::vector<Edge> path = {{0, 1}, {1, 2}};
  check(!(canonical_form(3, path, {1, 0, 0}).cert == canonical_form(3, path, {0, 1, 0}).cert),
        "colours enter the certificate");
}

void test_twins() {
  // Bowtie: the two free corners of each triangle are twins.
  const std::vector<Edge> bowtie = {{0, 1}, {0, 2}, {1, 2}, {2, 3}, {2, 4}, {3, 4}};
  check(twin_pairs(5, bowtie).size() == 2, "bowtie has two twin pairs");
}

void test_cluster_counts() {
  // Topologically distinct kagome clusters with 1..8 triangles (paper, p. 2).
  const Expansion kag = build_expansion("kagome", 8, 4);
  const std::vector<int> paper = {1, 1, 1, 1, 2, 2, 5, 7, 15};
  for (int n = 0; n <= 8; ++n) {
    check(kag.classes_per_order[n] == paper[n], "kagome classes at order " + std::to_string(n) + ": got " +
                                                    std::to_string(kag.classes_per_order[n]));
  }
  // The hexagon of six triangles has 12 spins; the largest 8-triangle trees 17.
  int hexagon = 0, max8 = 0;
  for (const auto& c : kag.classes) {
    if (c.order == 6 && c.nspins == 12) ++hexagon;
    if (c.order == 8) max8 = std::max(max8, c.nspins);
  }
  check(hexagon == 1, "one 12-site six-triangle cluster");
  check(max8 == 17, "largest eight-triangle cluster has 17 spins");

  // Site clusters of the square lattice per site are the fixed polyominoes.
  const Expansion sq = build_expansion("square", 9, 4);
  const std::vector<long long> polyominoes = {0, 1, 2, 6, 19, 63, 216, 760, 2725, 9910};
  for (int n = 1; n <= 9; ++n) check(sq.embeddings_per_cell[n] == polyominoes[n], "square fixed clusters");

  // The point groups give the same lattice constants as translations alone.
  for (const std::string lattice : {"kagome", "triangular", "square"}) {
    const int order = lattice == "square" ? 7 : 5;
    const Expansion a = build_expansion(lattice, order, 1);
    const Expansion b = build_expansion_translations_only(lattice, order);
    check(a.classes.size() == b.classes.size(), lattice + ": class count independent of point group");
    check(a.embeddings_per_cell == b.embeddings_per_cell, lattice + ": embeddings independent of point group");
    // Classes are discovered in different orders; their canonical structure
    // identifies them.
    std::map<std::pair<std::vector<Edge>, std::vector<std::vector<int>>>, double> constants;
    for (const auto& c : a.classes) constants[{c.bonds, c.blocks}] = c.lattice_constant;
    for (const auto& c : b.classes) {
      auto it = constants.find({c.bonds, c.blocks});
      check(it != constants.end() && close(it->second, c.lattice_constant, 1e-12),
            lattice + ": lattice constant independent of point group");
    }
  }
  // Every weight must vanish for the sum of all lattice constants weighted by
  // subcluster counts; a cheap consistency check is that each subcluster
  // precedes its parent.
  for (const auto* ex : {&kag, &sq}) {
    for (std::size_t c = 0; c < ex->classes.size(); ++c) {
      for (auto [s, count] : ex->classes[c].subclusters) {
        check(static_cast<std::size_t>(s) < c && count > 0, "subclusters precede parents");
      }
    }
  }
}

void compare_spectra(int n, const std::vector<Edge>& bonds, const Model& m, const std::string& what) {
  const Spectrum a = spectrum(n, bonds, m);
  const Spectrum b = spectrum_without_symmetries(n, bonds, m);
  std::vector<double> ea, eb;
  for (std::size_t k = 0; k < a.energy.size(); ++k) ea.insert(ea.end(), std::lround(a.degeneracy[k]), a.energy[k]);
  for (std::size_t k = 0; k < b.energy.size(); ++k) eb.insert(eb.end(), std::lround(b.degeneracy[k]), b.energy[k]);
  bool ok = ea.size() == eb.size();
  for (std::size_t k = 0; ok && k < ea.size(); ++k) ok = close(ea[k], eb[k], 1e-9);
  check(ok, "spectrum with symmetries matches brute force: " + what);
}

void test_ed() {
  const std::vector<Edge> bowtie = {{0, 1}, {0, 2}, {1, 2}, {2, 3}, {2, 4}, {3, 4}};
  // Chain of three triangles plus a dangling triangle (twins and odd n).
  const std::vector<Edge> tree = {{0, 1}, {0, 2}, {1, 2}, {2, 3}, {2, 4}, {3, 4}, {4, 5},
                                  {4, 6}, {5, 6}, {3, 7}, {3, 8}, {7, 8}};
  // Square plaquette with a twin pair of pendant sites (even n, flip parity).
  const std::vector<Edge> plaquette = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 4}, {0, 5}};
  const std::vector<Model> models = {{1, 0.0, 0.0, 0.0},   {1, 1.0, 0.0, 0.0}, {1, 0.37, 0.0, 0.0},
                                     {1, 0.37, 0.0, 0.21}, {1, 0.0, 0.5, 0.0}, {1, 0.0, 0.5, 0.25},
                                     {1, 0.3, 0.4, 0.0},   {1, 0.3, 0.4, 0.2}};
  for (std::size_t k = 0; k < models.size(); ++k) {
    compare_spectra(5, bowtie, models[k], "bowtie model " + std::to_string(k));
    compare_spectra(9, tree, models[k], "triangle tree model " + std::to_string(k));
    compare_spectra(6, plaquette, models[k], "plaquette model " + std::to_string(k));
  }
  // Heisenberg triangle: S = 3/2 at +3/4, two S = 1/2 doublets at -3/4.
  const Spectrum tri = spectrum(3, {{0, 1}, {1, 2}, {0, 2}}, {1, 1.0, 0.0, 0.0});
  check(tri.energy.size() == 2 && close(tri.energy[0], -0.75, 1e-12) && close(tri.degeneracy[0], 4, 1e-12) &&
            close(tri.energy[1], 0.75, 1e-12),
        "Heisenberg triangle spectrum");
}

void test_thermo() {
  // A free spin in a field h: Z = 2 cosh(h/2T).
  const Spectrum spin = spectrum(1, {}, {1, 0.0, 0.0, 0.8});
  const std::vector<double> temps = {0.05, 0.3, 1.0, 7.0};
  const Thermo th = thermodynamics(spin, temps);
  for (std::size_t t = 0; t < temps.size(); ++t) {
    const double x = 0.4 / temps[t];
    check(close(th.energy[t], -0.4 * std::tanh(x), 1e-12), "free spin energy");
    check(close(th.entropy[t], std::log(2 * std::cosh(x)) - x * std::tanh(x), 1e-12), "free spin entropy");
    check(close(th.heat[t], x * x / (std::cosh(x) * std::cosh(x)), 1e-12), "free spin heat");
  }
}

void test_resummation() {
  // ln 2 = 1 - 1/2 + 1/3 - ...
  std::vector<double> terms, partial;
  double sum = 0.0;
  for (int k = 0; k < 13; ++k) {
    terms.push_back((k % 2 ? -1.0 : 1.0) / (k + 1));
    sum += terms.back();
    partial.push_back(sum);
  }
  const double ln2 = std::log(2.0);
  check(std::abs(partial.back() - ln2) > 1e-2, "bare sum is slow");
  check(close(wynn(partial, 4), ln2, 1e-7), "Wynn accelerates ln 2");
  check(close(brezinski(partial, 3), ln2, 1e-6), "Brezinski accelerates ln 2");
  check(close(euler(terms, 2), ln2, 1e-4), "Euler accelerates ln 2");
  // Geometric partial sums are summed exactly by one Wynn cycle (Shanks).
  std::vector<double> geo = {1.0, 1.5, 1.75};
  check(close(wynn(geo, 1), 2.0, 1e-14), "Shanks on a geometric series");
  // Constant sequences are returned unchanged.
  std::vector<double> flat(7, 0.25);
  check(wynn(flat, 3) == 0.25 && brezinski(flat, 2) == 0.25, "converged sequences pass through");
}

void test_kagome_ising_entropy() {
  // Pauling: ln 2 + (2/3) ln(3/4) = 0.50136; the six-triangle hexagon brings
  // it to 0.50182 against the exact 0.50183 (paper, p. 2).
  const Expansion ex = build_expansion("kagome", 6, 4);
  const auto res = run_nlc(ex, {{1, 0.0, 0.0, 0.0}}, {1e-3}, {4, {}});
  double s5 = 0.0, s6 = 0.0;
  for (int n = 0; n <= 6; ++n) {
    if (n <= 5) s5 += res[0].entropy[n][0];
    s6 += res[0].entropy[n][0];
  }
  check(close(s5, std::log(2.0) + 2.0 / 3.0 * std::log(0.75), 1e-9), "Pauling entropy through five triangles");
  check(close(s6, 0.50182, 5e-6), "kagome Ising entropy through six triangles: " + std::to_string(s6));
  std::printf("kagome Ising S(T->0): five triangles %.5f, six triangles %.5f\n", s5, s6);

  // High temperature: S -> ln 2 and E -> 0 per site at every order.
  const auto hot = run_nlc(ex, {{1, 0.3, 0.2, 0.1}}, {1e6}, {4, {}});
  double s = 0.0;
  for (int n = 0; n <= 6; ++n) s += hot[0].entropy[n][0];
  check(close(s, std::log(2.0), 1e-9), "infinite-temperature entropy");
}

}  // namespace

int main() {
  use_single_threaded_lapack();
  test_canonical_form();
  test_twins();
  test_cluster_counts();
  test_ed();
  test_thermo();
  test_resummation();
  test_kagome_ising_entropy();
  if (failures) {
    std::fprintf(stderr, "%d failure(s)\n", failures);
    return 1;
  }
  std::printf("all tests passed\n");
  return 0;
}
