// Thermodynamics of a finite cluster from its spectrum.
#pragma once

#include "ed.hpp"

#include <vector>

namespace nlc {

// Extensive cluster totals on a temperature grid (k_B = 1).
struct Thermo {
  std::vector<double> energy;   // <H>
  std::vector<double> entropy;  // ln Z + <H>/T
  std::vector<double> heat;     // (<H^2> - <H>^2)/T^2
};

// Energies are measured from the ground state, so the entropy and specific
// heat stay accurate at temperatures far below the level spacing.
Thermo thermodynamics(const Spectrum& spec, const std::vector<double>& temps);

}  // namespace nlc
