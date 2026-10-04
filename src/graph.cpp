#include "graph.hpp"

#include <algorithm>
#include <functional>
#include <stdexcept>

namespace nlc {

namespace {

using Colors = std::array<int, kMaxVertices>;
using Adjacency = std::array<std::uint32_t, kMaxVertices>;

int pair_index(int n, int i, int j) { return i * n - i * (i + 1) / 2 + (j - i - 1); }

Adjacency adjacency(int n, const std::vector<Edge>& edges) {
  if (n > kMaxVertices) throw std::invalid_argument("graph too large for canonical labelling");
  Adjacency adj{};
  for (auto [u, v] : edges) {
    if (u == v || u < 0 || v < 0 || u >= n || v >= n) throw std::invalid_argument("bad edge");
    adj[u] |= 1u << v;
    adj[v] |= 1u << u;
  }
  return adj;
}

// Equitable refinement: a vertex's new colour is the rank of
// (colour, degree, sorted neighbour colours) among all vertices.
void refine(int n, const Adjacency& adj, Colors& colors) {
  struct Sig {
    std::array<int, kMaxVertices + 2> key;
    int len;
    int vertex;
  };
  std::array<Sig, kMaxVertices> sig;
  auto less = [](const Sig& x, const Sig& y) {
    return std::lexicographical_compare(x.key.begin(), x.key.begin() + x.len, y.key.begin(), y.key.begin() + y.len);
  };
  while (true) {
    for (int i = 0; i < n; ++i) {
      Sig& s = sig[i];
      s.vertex = i;
      s.key[0] = colors[i];
      int deg = 0;
      for (std::uint32_t m = adj[i]; m; m &= m - 1) s.key[2 + deg++] = colors[__builtin_ctz(m)];
      std::sort(s.key.begin() + 2, s.key.begin() + 2 + deg);
      s.key[1] = deg;
      s.len = deg + 2;
    }
    std::array<int, kMaxVertices> order;
    for (int i = 0; i < n; ++i) order[i] = i;
    std::sort(order.begin(), order.begin() + n, [&](int a, int b) { return less(sig[a], sig[b]); });
    Colors next{};
    int rank = 0;
    for (int k = 0; k < n; ++k) {
      if (k > 0 && less(sig[order[k - 1]], sig[order[k]])) ++rank;
      next[order[k]] = rank;
    }
    bool stable = true;
    for (int i = 0; i < n && stable; ++i) stable = next[i] == colors[i];
    colors = next;
    if (stable) return;
  }
}

struct Search {
  int n;
  Adjacency adj;
  bool have = false;
  std::array<std::uint64_t, kCertWords> best{};
  std::vector<int> best_perm;

  void leaf(const Colors& colors) {
    std::array<std::uint64_t, kCertWords> bits{};
    for (int i = 0; i < n; ++i) {
      for (std::uint32_t m = adj[i]; m; m &= m - 1) {
        const int a = colors[i], b = colors[__builtin_ctz(m)];
        if (a < b) {
          const int k = pair_index(n, a, b);
          bits[k / 64] |= std::uint64_t{1} << (63 - k % 64);
        }
      }
    }
    if (!have || bits < best) {
      have = true;
      best = bits;
      best_perm.assign(colors.begin(), colors.begin() + n);
    }
  }

  void recurse(Colors colors) {
    refine(n, adj, colors);
    std::array<int, kMaxVertices> count{};
    for (int i = 0; i < n; ++i) ++count[colors[i]];
    int target = -1;
    for (int c = 0; c < n && target < 0; ++c) {
      if (count[c] > 1) target = c;
    }
    if (target < 0) {
      leaf(colors);
      return;
    }
    for (int v = 0; v < n; ++v) {
      if (colors[v] != target) continue;
      Colors next = colors;
      for (int i = 0; i < n; ++i) {
        if (colors[i] > target || (colors[i] == target && i != v)) ++next[i];
      }
      recurse(next);
    }
  }
};

}  // namespace

std::size_t CertHash::operator()(const Cert& c) const {
  std::size_t h = static_cast<std::size_t>(c.n) * 0xc2b2ae3d27d4eb4fULL;
  for (std::uint64_t w : c.bits) h = (h ^ w) * 0x9e3779b97f4a7c15ULL + (h >> 29);
  for (int i = 0; i < c.n; ++i) h = (h ^ c.colors[i]) * 0x100000001b3ULL;
  return h;
}

Canon canonical_form(int n, const std::vector<Edge>& edges, const std::vector<int>& colors) {
  if (n == 0) return {};
  if (!colors.empty() && static_cast<int>(colors.size()) != n) throw std::invalid_argument("colour count mismatch");
  Search search{n, adjacency(n, edges), false, {}, {}};
  Colors initial{};
  for (int i = 0; i < static_cast<int>(colors.size()); ++i) {
    if (colors[i] < 0 || colors[i] > 255) throw std::invalid_argument("colour out of range");
    initial[i] = colors[i];
  }
  search.recurse(initial);

  Canon out;
  out.perm = search.best_perm;
  out.cert.n = n;
  out.cert.bits = search.best;
  for (int v = 0; v < n; ++v) out.cert.colors[out.perm[v]] = static_cast<std::uint8_t>(colors.empty() ? 0 : colors[v]);
  return out;
}

std::vector<Edge> twin_pairs(int n, const std::vector<Edge>& edges) {
  const Adjacency adj = adjacency(n, edges);
  std::uint32_t used = 0;
  std::vector<Edge> out;
  for (int u = 0; u < n; ++u) {
    if (used >> u & 1u) continue;
    for (int v = u + 1; v < n; ++v) {
      if (used >> v & 1u) continue;
      if ((adj[u] & ~(1u << v)) == (adj[v] & ~(1u << u))) {
        out.emplace_back(u, v);
        used |= (1u << u) | (1u << v);
        break;
      }
    }
  }
  return out;
}

}  // namespace nlc
