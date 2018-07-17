#include "surface_hopping.h"

#include <cassert>

////////////////////////////////////////////////////////////////////////////////

namespace pot {

//----------------------------------------------------------------------------//

double surface_hopping::force(size_t ndim, size_t natom, size_t nbead,
                              const double *crd, double *f) {
  assert(nbead > 0);
  assert(init_pot == true);
  assert(init_hop == true);

  size_t ndof = ndim * natom;
  assert(ndof > 0);

  // std::cerr << crd[0] << std::endl;

  check_allocation(nbead, state_id);

  double Eactive(0);
  double fAA[ndof];
  double fBB[ndof];
  double dF[ndof];
  for (size_t n = 0; n < nbead; ++n) {
    std::fill(fAA, fAA + ndof, 0.0);
    std::fill(fBB, fBB + ndof, 0.0);
    double EAA = VAA(crd + ndof * n, fAA);
    double EBB = VBB(crd + ndof * n, fBB);

    double dE = EBB - EAA;

    if (state_id[n] == 0) {
      std::copy(fAA, fAA + ndof, f + ndof * n);
      Eactive += EAA;
    } else {
      std::copy(fBB, fBB + ndof, f + ndof * n);
      Eactive += EBB;
    }

    // FIXME Not implemented
    // double EBath = bath_force(crd + ndof*n, f + ndof*n);
    // Eactive += EBath;

    // Determine if the active state of this bead should be changed
    double random_number = ((double)rand() / (RAND_MAX));
    assert(random_number >= 0.0 && random_number <= 1.0);

    if (hop_probability(dE, state_id[n]) > random_number) {
      // Hop successful, change state id of this bead
      if (state_id[n] == 0)
        state_id[n] = 1;
      else
        state_id[n] = 0;
    }
  }

  return Eactive;
}

//----------------------------------------------------------------------------//

double surface_hopping::fermi_function(const double dE, const double mu) {
  double exp_arg = beta * (dE - mu);
  if (exp_arg > 100.0)
    return 0.0;
  else
    return 1.0 / (1.0 + std::exp(exp_arg));
}

//----------------------------------------------------------------------------//

double surface_hopping::hop_probability(const double dE, int active_state) {
  double fE_L = fermi_function(dE, voltage / 2.0);
  double fE_R = fermi_function(dE, -voltage / 2.0);

  double Gamma_L = Gamma / 2.0;
  double Gamma_R = Gamma / 2.0;

  double f_bar = (Gamma_L * fE_L + Gamma_R * fE_R) / Gamma;

  if (active_state == 0)
    return Gamma * dt * f_bar;
  else
    return Gamma * dt * (1.0 - f_bar);
}

//----------------------------------------------------------------------------//

// These parameters define the hopping probability
void surface_hopping::set_hopping_params(double *params) {
  Gamma = params[0];
  dt = params[1];
  beta = params[2];
  voltage = params[3];

  // prng.seed(19104);

  init_hop = true;
}

//----------------------------------------------------------------------------//

void surface_hopping::print_state_params() {
  std::cout << "# surface_hopping_dynamics" << std::endl;
  std::cout << "# GammaEl = " << Gamma << std::endl;
  std::cout << "# V     = " << voltage << std::endl;
}

//----------------------------------------------------------------------------//

} // namespace pot

////////////////////////////////////////////////////////////////////////////////
