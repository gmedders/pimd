#include <cassert>
#include <cmath>

#include "nhc.h"

namespace parts {
namespace nhc {

//
// thermostat is stored as eta1 ... etaM, vel1 ... velM, acc1 ... accM
//

size_t size(size_t M) { return 3 * M; }

void initialize(size_t M, double *thermo, const double &tau,
                std::mt19937 &prng) {
  assert(M > 1);
  assert(thermo != 0);
  assert(tau > 0.0);

  double *eta = thermo;
  double *vel = eta + M;
  double *acc = vel + M;

  const double tau1 = 1.0 / tau; // == sqrt(kT/Q), Q = kT*tau*tau

  std::normal_distribution<> rand_01{0, 1};

  for (size_t i = 0; i < M; ++i) {
    eta[i] = 0.0;
    vel[i] = tau1 * rand_01(prng);
    if (i > 0)
      acc[i] = vel[i - 1] * vel[i - 1] - tau1 * tau1;
    // acc[0] is computed in advance()
  }
}

//
// follows eq (4.9)
//

double advance(size_t M, double *thermo, const double &tau, const double &Ek2kT,
               const double &dt) {
  assert(M > 1 && thermo != 0 && tau > 0.0);

  double *eta = thermo;
  double *vel = eta + M;
  double *acc = vel + M;

  const double dt2 = dt / 2;
  const double dt4 = dt / 4;

  const double tau2 = 1.0 / tau / tau;

  // force acting on the first oscillator
  acc[0] = tau2 * (Ek2kT - 1.0);

  // update the velocities/positions
  size_t m = M - 1;
  vel[m] += dt2 * acc[m];
  eta[m] += dt * vel[m];

  while (m != 0) {
    const double tmp = std::exp(-dt4 * vel[m--]);
    vel[m] = tmp * (tmp * vel[m] + dt2 * acc[m]);
    eta[m] += dt * vel[m];
  }

  const double aa = std::exp(-dt * vel[0]); // velocity scale factor

  // update velocities
  acc[0] = tau2 * (Ek2kT * aa * aa - 1.0);
  const double tmp = std::exp(-dt4 * vel[1]);
  vel[0] = tmp * (tmp * vel[0] + dt2 * acc[0]);

  m = 0;
  while (m + 2 < M) {
    acc[m + 1] = vel[m] * vel[m] - tau2;
    ++m;
    const double tmp = std::exp(-dt4 * vel[m + 1]);
    vel[m] = tmp * (tmp * vel[m] + dt2 * acc[m]);
  }

  acc[m + 1] = vel[m] * vel[m] - tau2;
  ++m;
  vel[m] += dt2 * acc[m];

  return aa;
}

//
// returns sum of tau^2*vel^2/2 + eta
//

double invariant(size_t M, const double *thermo, const double &tau) {
  assert(M > 1 && thermo != 0);

  const double *eta = thermo;
  const double *vel = eta + M;

  double sum(0);
  while (M-- != 0) {
    const double tmp = tau * (*vel++);
    sum += tmp * tmp / 2 + *eta++;
  }

  return sum;
}

} // namespace nhc
} // namespace parts
