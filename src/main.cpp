// nlc: numerical linked-cluster expansion for spin-1/2 XXZ models in
// transverse and longitudinal fields (Rigol, Bryant and Singh, PRL 97,
// 187202 (2006)).
//
//   nlc clusters --lattice L --order N
//   nlc run --lattice L --order N --out PREFIX [--jz 1] [--jperp 0,0.25] [--hx 0] [--hz 0]
//           [--tmin 1e-3] [--tmax 100] [--nt 401] [--threads T] [--memory-gb 8] [--euler-start 2]
//
// L is kagome or triangular (triangle expansions, N triangles) or square
// (site expansion, N sites).  --jperp, --hx and --hz take comma-separated
// lists; `run` writes PREFIX_jp<J>_hx<hx>_hz<hz>.csv for every combination,
// with per-site energy E, entropy S and specific heat C on a logarithmic
// temperature grid: the bare partial sums through every order (X_n), Wynn
// with c cycles (X_wynn<c>), Brezinski with c cycles (X_brez<c>) and the
// Euler transformation (X_euler).
#include "expansion.hpp"
#include "nlc.hpp"
#include "resum.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

using namespace nlc;

namespace {

struct Args {
  std::string command;
  std::map<std::string, std::string> opts;

  bool has(const std::string& k) const { return opts.count(k) > 0; }
  std::string get(const std::string& k, const std::string& def = "") const {
    auto it = opts.find(k);
    return it == opts.end() ? def : it->second;
  }
  int integer(const std::string& k, int def) const { return has(k) ? std::stoi(get(k)) : def; }
  double real(const std::string& k, double def) const { return has(k) ? std::stod(get(k)) : def; }
  std::vector<double> reals(const std::string& k, double def) const {
    if (!has(k)) return {def};
    std::vector<double> out;
    std::stringstream ss(get(k));
    for (std::string item; std::getline(ss, item, ',');) out.push_back(std::stod(item));
    return out;
  }
};

Args parse(int argc, char** argv) {
  if (argc < 2) throw std::invalid_argument("usage: nlc {clusters|run} [options]");
  Args a;
  a.command = argv[1];
  for (int i = 2; i < argc; ++i) {
    std::string k = argv[i];
    if (k.rfind("--", 0) != 0 || i + 1 >= argc) throw std::invalid_argument("bad option " + k);
    a.opts[k.substr(2)] = argv[++i];
  }
  return a;
}

std::string fmt(double x) {
  char buf[32];
  std::snprintf(buf, sizeof buf, "%g", x);
  return buf;
}

int first_order(const Expansion& ex) { return ex.classes.front().order; }

void print_clusters(const Expansion& ex) {
  std::printf("%s lattice, orders %d..%d, %zu classes\n", ex.lattice.c_str(), first_order(ex), ex.max_order,
              ex.classes.size());
  std::printf("%6s %8s %14s %10s\n", "order", "classes", "per unit cell", "max spins");
  for (int n = first_order(ex); n <= ex.max_order; ++n) {
    int spins = 0;
    for (const auto& c : ex.classes) {
      if (c.order == n) spins = std::max(spins, c.nspins);
    }
    std::printf("%6d %8d %14lld %10d\n", n, ex.classes_per_order[n], ex.embeddings_per_cell[n], spins);
  }
}

void write_csv(const std::string& path, const Expansion& ex, const NlcResult& res, int euler_start) {
  const int n0 = first_order(ex), nmax = ex.max_order;
  const int len = nmax - n0 + 1;
  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot write " + path);
  out.precision(12);
  const std::vector<std::pair<std::string, const std::vector<std::vector<double>>*>> props = {
      {"E", &res.energy}, {"S", &res.entropy}, {"C", &res.heat}};
  out << "T";
  for (const auto& [name, _] : props) {
    for (int n = n0; n <= nmax; ++n) out << ',' << name << '_' << n;
    for (int c = 1; 2 * c + 1 <= len; ++c) out << ',' << name << "_wynn" << c;
    for (int c = 1; 3 * c + 1 <= len; ++c) out << ',' << name << "_brez" << c;
    if (euler_start - n0 < len) out << ',' << name << "_euler";
  }
  out << '\n';
  for (std::size_t t = 0; t < res.temps.size(); ++t) {
    out << res.temps[t];
    for (const auto& [name, terms] : props) {
      std::vector<double> term(len), partial(len);
      double sum = 0.0;
      for (int n = n0; n <= nmax; ++n) {
        term[n - n0] = (*terms)[n][t];
        sum += term[n - n0];
        partial[n - n0] = sum;
        out << ',' << sum;
      }
      for (int c = 1; 2 * c + 1 <= len; ++c) out << ',' << wynn(partial, c);
      for (int c = 1; 3 * c + 1 <= len; ++c) out << ',' << brezinski(partial, c);
      if (euler_start - n0 < len) out << ',' << euler(term, std::max(euler_start - n0, 0));
    }
    out << '\n';
  }
}

int run(const Args& args) {
  const std::string lattice = args.get("lattice", "kagome");
  const int order = args.integer("order", 6);
  const int threads = args.integer("threads", static_cast<int>(std::max(1u, std::thread::hardware_concurrency())));
  auto log = [](const std::string& s) { std::cerr << s << std::endl; };

  auto t0 = std::chrono::steady_clock::now();
  const Expansion ex = build_expansion(lattice, order, threads);
  const double tgeo = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

  if (args.command == "clusters") {
    print_clusters(ex);
    std::printf("built in %.2f s\n", tgeo);
    return 0;
  }
  if (args.command != "run") throw std::invalid_argument("unknown command " + args.command);
  if (!args.has("out")) throw std::invalid_argument("run needs --out PREFIX");

  std::vector<Model> models;
  const double jz = args.real("jz", 1.0);
  for (double jp : args.reals("jperp", 0.0)) {
    for (double hx : args.reals("hx", 0.0)) {
      for (double hz : args.reals("hz", 0.0)) models.push_back({jz, jp, hx, hz});
    }
  }
  const double tmin = args.real("tmin", 1e-3), tmax = args.real("tmax", 100.0);
  const int nt = args.integer("nt", 401);
  std::vector<double> temps(nt);
  for (int i = 0; i < nt; ++i) temps[i] = tmin * std::pow(tmax / tmin, nt == 1 ? 0.0 : static_cast<double>(i) / (nt - 1));

  set_memory_budget(static_cast<std::size_t>(args.real("memory-gb", 8.0) * (1ull << 30)));
  use_single_threaded_lapack();
  log("expansion: " + lattice + " to order " + std::to_string(order) + ", " + std::to_string(ex.classes.size()) +
      " classes (" + fmt(tgeo) + " s); " + std::to_string(models.size()) + " models, " + std::to_string(threads) +
      " threads");
  t0 = std::chrono::steady_clock::now();
  const auto results = run_nlc(ex, models, temps, {threads, log});
  log("diagonalisation: " + fmt(std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count()) + " s");

  const int euler_start = args.integer("euler-start", 2);
  for (const auto& res : results) {
    const std::string path = args.get("out") + "_jp" + fmt(res.model.jperp) + "_hx" + fmt(res.model.hx) + "_hz" +
                             fmt(res.model.hz) + ".csv";
    write_csv(path, ex, res, euler_start);
    log("wrote " + path + " (largest block " + std::to_string(res.largest_block) + ")");
  }
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    return run(parse(argc, argv));
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << '\n';
    return 1;
  }
}
