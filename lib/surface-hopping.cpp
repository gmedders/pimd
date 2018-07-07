#include "surface-hopping.h"

#include <cassert>

////////////////////////////////////////////////////////////////////////////////

namespace pot {

//----------------------------------------------------------------------------//

// Both determines whether the vector is the correct size
// and reallocates if necessary
void surface_hopping::check_allocation(size_t n_elements, arma::ivec &myvec) {
  if (myvec.n_elem != n_elements)
    myvec.set_size(n_elements);
}

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

void surface_hopping::set_all_bead_states(int initial_state_id, int nbead) {
  assert(initial_state_id == 0 || initial_state_id == 1);

  m_nbead = nbead;
  // Set for all beads
  check_allocation(m_nbead, state_id);
  for (int n = 0; n < m_nbead; ++n) {
    state_id[n] = initial_state_id;
  }
}

//----------------------------------------------------------------------------//

void surface_hopping::set_individual_bead_states(
    std::vector<int> &initial_state_ids) {
  m_nbead = initial_state_ids.size();

  // Set for all beads
  check_allocation(m_nbead, state_id);
  for (int n = 0; n < m_nbead; ++n) {
    int id = initial_state_ids[n];
    assert(id == 0 || id == 1);
    state_id[n] = id;
  }
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

double surface_hopping::avg_active_state() {
  return sum_active_state() / ((double)m_nbead);
}

//----------------------------------------------------------------------------//

double surface_hopping::sum_active_state() {
  double sum(0);
  for (int i = 0; i < m_nbead; ++i)
    sum += state_id[i];
  return sum;
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
