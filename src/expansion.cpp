#include "expansion.hpp"

#include "lattice.hpp"

#include <boost/unordered/unordered_flat_map.hpp>

#include <algorithm>
#include <array>
#include <atomic>
#include <map>
#include <stdexcept>
#include <thread>

namespace nlc {

namespace {

using ClassMap = boost::unordered_flat_map<Cert, int, CertHash>;

enum class Family { triangles, sites };

struct Spec {
  Family family;
  Lattice blocks;     // lattice of building blocks (sites for the site expansion)
  int spins_per_cell;  // spins of the physical lattice per unit cell of `blocks`
  std::array<Site, 3> (*corners)(Site);
};

std::array<Site, 3> kagome_corners(Site s) {
  // Kagome spins are the midpoints of honeycomb bonds, keyed by s + t; the
  // triangle at honeycomb site s holds its three bonds.
  const auto nb = honeycomb_lattice().neighbours(s);
  if (nb.size() != 3) throw std::logic_error("honeycomb site without three neighbours");
  return {s + nb[0], s + nb[1], s + nb[2]};
}

std::array<Site, 3> up_triangle_corners(Site s) { return {s, s + Site{1, 0}, s + Site{0, 1}}; }

Spec spec_for(const std::string& lattice) {
  if (lattice == "kagome") return {Family::triangles, honeycomb_lattice(), 3, kagome_corners};
  if (lattice == "triangular") return {Family::triangles, up_triangle_lattice(), 1, up_triangle_corners};
  if (lattice == "square") return {Family::sites, square_lattice(), 1, nullptr};
  throw std::invalid_argument("unknown lattice " + lattice + " (kagome, triangular, square)");
}

// Incidence graph of a block cluster: spins 0..ns-1 (colour 0) followed by
// blocks (colour 1).
Canon block_canon(int nspins, const std::vector<std::vector<int>>& blocks) {
  std::vector<Edge> edges;
  std::vector<int> colors(nspins, 0);
  for (std::size_t b = 0; b < blocks.size(); ++b) {
    colors.push_back(1);
    for (int v : blocks[b]) edges.emplace_back(v, nspins + static_cast<int>(b));
  }
  return canonical_form(nspins + static_cast<int>(blocks.size()), edges, colors);
}

std::vector<Edge> block_bonds(const std::vector<std::vector<int>>& blocks) {
  std::vector<Edge> bonds;
  for (const auto& blk : blocks) {
    for (std::size_t i = 0; i < blk.size(); ++i) {
      for (std::size_t j = i + 1; j < blk.size(); ++j) bonds.emplace_back(std::min(blk[i], blk[j]), std::max(blk[i], blk[j]));
    }
  }
  std::sort(bonds.begin(), bonds.end());
  return bonds;
}

// Relabels a block cluster (spins 0..ns-1) canonically.
ClusterClass make_block_class(int order, int nspins, const std::vector<std::vector<int>>& blocks, const Canon& canon) {
  ClusterClass cls;
  cls.order = order;
  cls.nspins = nspins;
  for (const auto& blk : blocks) {
    std::vector<int> b;
    for (int v : blk) b.push_back(canon.perm[v]);
    std::sort(b.begin(), b.end());
    cls.blocks.push_back(std::move(b));
  }
  std::sort(cls.blocks.begin(), cls.blocks.end());
  cls.bonds = block_bonds(cls.blocks);
  return cls;
}

ClusterClass make_site_class(int nspins, const std::vector<Edge>& edges, const Canon& canon) {
  ClusterClass cls;
  cls.order = nspins;
  cls.nspins = nspins;
  for (auto [u, v] : edges) {
    const int a = canon.perm[u], b = canon.perm[v];
    cls.bonds.emplace_back(std::min(a, b), std::max(a, b));
  }
  std::sort(cls.bonds.begin(), cls.bonds.end());
  return cls;
}

bool connected(std::uint32_t mask, const std::vector<std::uint32_t>& adj) {
  std::uint32_t seen = mask & (~mask + 1), frontier = seen;
  while (frontier) {
    std::uint32_t next = 0;
    for (std::uint32_t m = frontier; m; m &= m - 1) next |= adj[__builtin_ctz(m)];
    next &= mask & ~seen;
    seen |= next;
    frontier = next;
  }
  return seen == mask;
}

int lookup(const ClassMap& map, const Cert& cert) {
  auto it = map.find(cert);
  if (it == map.end()) throw std::logic_error("subcluster class missing from expansion");
  return it->second;
}

std::vector<std::pair<int, int>> block_subclusters(const ClusterClass& cls, const ClassMap& map) {
  const int m = static_cast<int>(cls.blocks.size());
  std::vector<std::uint32_t> adj(m, 0);
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < m; ++j) {
      if (i == j) continue;
      for (int v : cls.blocks[i]) {
        if (std::find(cls.blocks[j].begin(), cls.blocks[j].end(), v) != cls.blocks[j].end()) adj[i] |= 1u << j;
      }
    }
  }
  std::map<int, int> counts;
  counts[0] = cls.nspins;  // single sites
  for (std::uint32_t mask = 1; mask + 1 < (1u << m); ++mask) {
    if (!connected(mask, adj)) continue;
    std::vector<int> local(cls.nspins, -1);
    int ns = 0;
    std::vector<std::vector<int>> blocks;
    for (int b = 0; b < m; ++b) {
      if (!(mask >> b & 1u)) continue;
      std::vector<int> blk;
      for (int v : cls.blocks[b]) {
        if (local[v] < 0) local[v] = ns++;
        blk.push_back(local[v]);
      }
      blocks.push_back(std::move(blk));
    }
    ++counts[lookup(map, block_canon(ns, blocks).cert)];
  }
  return {counts.begin(), counts.end()};
}

std::vector<std::pair<int, int>> site_subclusters(const ClusterClass& cls, const ClassMap& map) {
  const int n = cls.nspins;
  std::vector<std::uint32_t> adj(n, 0);
  for (auto [i, j] : cls.bonds) {
    adj[i] |= 1u << j;
    adj[j] |= 1u << i;
  }
  std::map<int, int> counts;
  for (std::uint32_t mask = 1; mask + 1 < (1u << n); ++mask) {
    if (!connected(mask, adj)) continue;
    std::vector<int> local(n, -1);
    int ns = 0;
    for (int v = 0; v < n; ++v) {
      if (mask >> v & 1u) local[v] = ns++;
    }
    std::vector<Edge> edges;
    for (auto [i, j] : cls.bonds) {
      if (local[i] >= 0 && local[j] >= 0) edges.emplace_back(local[i], local[j]);
    }
    ++counts[lookup(map, canonical_form(ns, edges).cert)];
  }
  return {counts.begin(), counts.end()};
}

Expansion build(const std::string& lattice, const Spec& spec, int max_order, int threads) {
  if (max_order < 0) throw std::invalid_argument("negative order");
  Expansion ex;
  ex.lattice = lattice;
  ex.max_order = max_order;
  ex.classes_per_order.assign(max_order + 1, 0);
  ex.embeddings_per_cell.assign(max_order + 1, 0);
  ClassMap map;
  const double per_site = 1.0 / spec.spins_per_cell;

  if (spec.family == Family::triangles) {
    const Canon site = block_canon(1, {});
    map.emplace(site.cert, 0);
    ClusterClass single;
    single.nspins = 1;
    single.lattice_constant = 1.0;
    ex.classes.push_back(single);
    ex.classes_per_order[0] = 1;
    ex.embeddings_per_cell[0] = spec.spins_per_cell;
  }
  if (max_order >= 1) {
    const auto reps = enumerate_site_clusters(spec.blocks, max_order);
    for (int order = 1; order <= max_order; ++order) {
      for (const auto& cluster : reps[order]) {
        const int orbit = site_orbit_size(spec.blocks, cluster);
        ex.embeddings_per_cell[order] += orbit;
        Cert cert;
        ClusterClass fresh;
        if (spec.family == Family::triangles) {
          std::vector<Site> spins;
          for (const auto& b : cluster) {
            for (const auto& c : spec.corners(b)) spins.push_back(c);
          }
          std::sort(spins.begin(), spins.end());
          spins.erase(std::unique(spins.begin(), spins.end()), spins.end());
          std::vector<std::vector<int>> blocks;
          for (const auto& b : cluster) {
            std::vector<int> blk;
            for (const auto& c : spec.corners(b)) {
              blk.push_back(static_cast<int>(std::lower_bound(spins.begin(), spins.end(), c) - spins.begin()));
            }
            blocks.push_back(std::move(blk));
          }
          const int ns = static_cast<int>(spins.size());
          const Canon canon = block_canon(ns, blocks);
          cert = canon.cert;
          if (!map.contains(cert)) fresh = make_block_class(order, ns, blocks, canon);
        } else {
          const int ns = static_cast<int>(cluster.size());
          std::vector<Edge> edges;
          for (const auto& [p, q] : induced_bonds(spec.blocks, cluster)) {
            edges.emplace_back(std::lower_bound(cluster.begin(), cluster.end(), p) - cluster.begin(),
                               std::lower_bound(cluster.begin(), cluster.end(), q) - cluster.begin());
          }
          const Canon canon = canonical_form(ns, edges);
          cert = canon.cert;
          if (!map.contains(cert)) fresh = make_site_class(ns, edges, canon);
        }
        auto [it, inserted] = map.try_emplace(cert, static_cast<int>(ex.classes.size()));
        if (inserted) {
          ex.classes.push_back(std::move(fresh));
          ++ex.classes_per_order[order];
        }
        ex.classes[it->second].lattice_constant += orbit * per_site;
      }
    }
  }

  std::atomic<std::size_t> next{0};
  auto work = [&] {
    for (std::size_t c; (c = next++) < ex.classes.size();) {
      auto& cls = ex.classes[c];
      if (cls.nspins == 1) continue;
      cls.subclusters = spec.family == Family::triangles ? block_subclusters(cls, map) : site_subclusters(cls, map);
    }
  };
  {
    std::vector<std::jthread> pool;
    for (int t = 1; t < std::max(threads, 1); ++t) pool.emplace_back(work);
    work();
  }
  return ex;
}

}  // namespace

Expansion build_expansion(const std::string& lattice, int max_order, int threads) {
  return build(lattice, spec_for(lattice), max_order, threads);
}

Expansion build_expansion_translations_only(const std::string& lattice, int max_order) {
  Spec spec = spec_for(lattice);
  spec.blocks = without_point_group(spec.blocks);
  return build(lattice, spec, max_order, 1);
}

}  // namespace nlc
