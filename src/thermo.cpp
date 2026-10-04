#include "thermo.hpp"

#include <cmath>
#include <stdexcept>

namespace nlc {

Thermo thermodynamics(const Spectrum& spec, const std::vector<double>& temps) {
  if (spec.energy.empty()) throw std::invalid_argument("empty spectrum");
  const std::size_t nt = temps.size(), nl = spec.energy.size();
  const double e0 = spec.energy.front();
  Thermo out;
  out.energy.resize(nt);
  out.entropy.resize(nt);
  out.heat.resize(nt);
  std::vector<double> boltz(nl);
  for (std::size_t t = 0; t < nt; ++t) {
    const double beta = 1.0 / temps[t];
    double z = 0.0, e1 = 0.0;
    for (std::size_t k = 0; k < nl; ++k) {
      boltz[k] = spec.degeneracy[k] * std::exp(-beta * (spec.energy[k] - e0));
      z += boltz[k];
      e1 += boltz[k] * (spec.energy[k] - e0);
    }
    e1 /= z;
    double var = 0.0;
    for (std::size_t k = 0; k < nl; ++k) {
      const double d = spec.energy[k] - e0 - e1;
      var += boltz[k] * d * d;
    }
    var /= z;
    out.energy[t] = e0 + e1;
    out.entropy[t] = std::log(z) + beta * e1;
    out.heat[t] = beta * beta * var;
  }
  return out;
}

}  // namespace nlc
