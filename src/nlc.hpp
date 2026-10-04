// Assembly of the linked-cluster sum from per-cluster properties.
#pragma once

#include "expansion.hpp"
#include "ed.hpp"

#include <functional>
#include <string>
#include <vector>

namespace nlc {

// values[c][t] is an extensive property of class c at grid point t.  Returns
// terms[n][t] = sum over classes c of order n of L(c) W(c) at t, with the
// weights W from the inclusion-exclusion recursion (Eq. 2).
std::vector<std::vector<double>> order_terms(const Expansion& ex, const std::vector<std::vector<double>>& values);

// Per-site energy, entropy and specific heat, resolved by order.
struct NlcResult {
  Model model;
  std::vector<double> temps;
  std::vector<std::vector<double>> energy, entropy, heat;  // [order][t]
  std::size_t largest_block = 0;
};

struct RunOptions {
  int threads = 1;
  std::function<void(const std::string&)> log;
};

// Diagonalises every class for every model (in parallel, largest clusters
// first) and assembles the sums.
std::vector<NlcResult> run_nlc(const Expansion& ex, const std::vector<Model>& models, const std::vector<double>& temps,
                               const RunOptions& opts);

}  // namespace nlc
