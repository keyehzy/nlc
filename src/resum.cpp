#include "resum.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace nlc {

namespace {

// True when the window has no variation beyond rounding.
bool converged(std::span<const double> window) {
  double scale = 0.0;
  for (double x : window) scale = std::max(scale, std::abs(x));
  const double last = window.back();
  for (double x : window) {
    if (std::abs(x - last) > 1e-13 * scale) return false;
  }
  return true;
}

double finite_or(double x, double fallback) { return std::isfinite(x) ? x : fallback; }

}  // namespace

double wynn(std::span<const double> partial, int cycles) {
  const std::size_t need = 2 * static_cast<std::size_t>(cycles) + 1;
  if (cycles < 1 || partial.size() < need) throw std::invalid_argument("too few partial sums for Wynn");
  const auto window = partial.subspan(partial.size() - need);
  if (converged(window)) return window.back();
  std::vector<double> prev(need + 1, 0.0), cur(window.begin(), window.end());
  for (std::size_t k = 1; k < need; ++k) {
    std::vector<double> next(cur.size() - 1);
    for (std::size_t j = 0; j + 1 < cur.size(); ++j) next[j] = prev[j + 1] + 1.0 / (cur[j + 1] - cur[j]);
    prev = std::move(cur);
    cur = std::move(next);
  }
  return finite_or(cur[0], window.back());
}

double brezinski(std::span<const double> partial, int cycles) {
  const std::size_t need = 3 * static_cast<std::size_t>(cycles) + 1;
  if (cycles < 1 || partial.size() < need) throw std::invalid_argument("too few partial sums for Brezinski");
  const auto window = partial.subspan(partial.size() - need);
  if (converged(window)) return window.back();
  std::vector<double> odd(need + 1, 0.0), even(window.begin(), window.end());
  for (int c = 0; c < cycles; ++c) {
    std::vector<double> next_odd(even.size() - 1);
    for (std::size_t j = 0; j < next_odd.size(); ++j) next_odd[j] = odd[j + 1] + 1.0 / (even[j + 1] - even[j]);
    std::vector<double> next_even(even.size() - 3);
    for (std::size_t j = 0; j < next_even.size(); ++j) {
      const double num = (even[j + 2] - even[j + 1]) * (next_odd[j + 2] - next_odd[j + 1]);
      const double den = next_odd[j + 2] - 2.0 * next_odd[j + 1] + next_odd[j];
      next_even[j] = even[j + 1] + num / den;
    }
    odd = std::move(next_odd);
    even = std::move(next_even);
  }
  return finite_or(even[0], window.back());
}

double euler(std::span<const double> terms, int start) {
  if (start < 0 || static_cast<std::size_t>(start) >= terms.size()) throw std::invalid_argument("bad Euler start");
  double sum = 0.0;
  for (int n = 0; n < start; ++n) sum += terms[n];
  std::vector<double> u(terms.begin() + start, terms.end());
  for (std::size_t k = 1; k < u.size(); k += 2) u[k] = -u[k];
  double scale = 0.5, sign = 1.0;
  while (!u.empty()) {
    sum += sign * scale * u[0];
    for (std::size_t k = 0; k + 1 < u.size(); ++k) u[k] = u[k + 1] - u[k];
    u.pop_back();
    scale *= 0.5;
    sign = -sign;
  }
  return sum;
}

}  // namespace nlc
