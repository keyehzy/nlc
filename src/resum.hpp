// Sequence extrapolation of the NLC partial sums (Guttmann; Numerical
// Recipes 5.1).  Every routine returns the last partial sum when the
// sequence has already converged to rounding or the table breaks down.
#pragma once

#include <span>

namespace nlc {

// Wynn's epsilon algorithm: epsilon_{2c} built from the last 2c + 1 partial
// sums, i.e. c cycles of improvement.
double wynn(std::span<const double> partial, int cycles);

// Brezinski's theta algorithm: theta_{2c} from the last 3c + 1 partial sums.
double brezinski(std::span<const double> partial, int cycles);

// Euler transformation: terms[0..start-1] summed directly, the remaining
// terms treated as an alternating series sum_k (-1)^k u_k with
// u_k = (-1)^k terms[start + k] and resummed as
// sum_s (-1)^s Delta^s u_0 / 2^(s+1).
double euler(std::span<const double> terms, int start);

}  // namespace nlc
