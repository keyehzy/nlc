#include "nlc.hpp"

#include "thermo.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <mutex>
#include <thread>

namespace nlc {

std::vector<std::vector<double>> order_terms(const Expansion& ex, const std::vector<std::vector<double>>& values) {
  const std::size_t nc = ex.classes.size();
  const std::size_t nt = values.empty() ? 0 : values[0].size();
  std::vector<std::vector<double>> weight(nc);
  std::vector<std::vector<double>> terms(ex.max_order + 1, std::vector<double>(nt, 0.0));
  for (std::size_t c = 0; c < nc; ++c) {
    const auto& cls = ex.classes[c];
    weight[c] = values[c];
    for (auto [s, count] : cls.subclusters) {
      for (std::size_t t = 0; t < nt; ++t) weight[c][t] -= count * weight[s][t];
    }
    for (std::size_t t = 0; t < nt; ++t) terms[cls.order][t] += cls.lattice_constant * weight[c][t];
  }
  return terms;
}

std::vector<NlcResult> run_nlc(const Expansion& ex, const std::vector<Model>& models, const std::vector<double>& temps,
                               const RunOptions& opts) {
  struct Task {
    int cls;
    int model;
  };
  std::vector<Task> tasks;
  for (int m = 0; m < static_cast<int>(models.size()); ++m) {
    for (int c = 0; c < static_cast<int>(ex.classes.size()); ++c) tasks.push_back({c, m});
  }
  auto quantum = [&](const Model& m) { return m.jperp != 0.0 || m.hx != 0.0; };
  // Largest Hilbert spaces first, and models without Sz conservation ahead of
  // those with it, so the long diagonalisations start early.
  std::stable_sort(tasks.begin(), tasks.end(), [&](const Task& a, const Task& b) {
    const auto key = [&](const Task& t) {
      const Model& m = models[t.model];
      return std::pair{quantum(m) ? ex.classes[t.cls].nspins : 0, m.hx != 0.0 ? 1 : 0};
    };
    return key(a) > key(b);
  });

  std::vector<std::vector<Thermo>> thermo(models.size(), std::vector<Thermo>(ex.classes.size()));
  std::vector<std::size_t> largest(models.size(), 0);
  std::mutex mutex;
  std::atomic<std::size_t> next{0};
  const auto start = std::chrono::steady_clock::now();
  auto work = [&] {
    for (std::size_t k; (k = next++) < tasks.size();) {
      const auto [c, m] = tasks[k];
      const auto& cls = ex.classes[c];
      const auto t0 = std::chrono::steady_clock::now();
      const Spectrum spec = spectrum(cls.nspins, cls.bonds, models[m]);
      thermo[m][c] = thermodynamics(spec, temps);
      const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      std::lock_guard lock(mutex);
      largest[m] = std::max(largest[m], spec.largest_block);
      if (opts.log && secs > 5.0) {
        const double total = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        char buf[160];
        std::snprintf(buf, sizeof buf, "  class %d (order %d, %d spins, model %d): block %zu, %.1f s [%.0f s, %zu/%zu]",
                      c, cls.order, cls.nspins, m, spec.largest_block, secs, total, k + 1, tasks.size());
        opts.log(buf);
      }
    }
  };
  {
    std::vector<std::jthread> pool;
    for (int t = 1; t < std::max(opts.threads, 1); ++t) pool.emplace_back(work);
    work();
  }

  std::vector<NlcResult> out;
  for (std::size_t m = 0; m < models.size(); ++m) {
    NlcResult res;
    res.model = models[m];
    res.temps = temps;
    res.largest_block = largest[m];
    std::vector<std::vector<double>> e, s, h;
    for (auto& th : thermo[m]) {
      e.push_back(std::move(th.energy));
      s.push_back(std::move(th.entropy));
      h.push_back(std::move(th.heat));
    }
    res.energy = order_terms(ex, e);
    res.entropy = order_terms(ex, s);
    res.heat = order_terms(ex, h);
    out.push_back(std::move(res));
  }
  return out;
}

}  // namespace nlc
