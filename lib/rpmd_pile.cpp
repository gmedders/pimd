#include <cassert>
#include <cmath>

#include <algorithm>

#include "rpmd_pile.h"

////////////////////////////////////////////////////////////////////////////////

namespace parts {

//----------------------------------------------------------------------------//

double rpmd_pile::calc_KE(void) {
  double e(0);
  for (size_t n = 0; n < nbeads(); ++n) {
    for (size_t i = 0; i < ndofs(); ++i) {
      double mass = m_mass(i);
      const double Ekin2 = m_mom_cart(i, n) * m_mom_cart(i, n) / mass;
      e += Ekin2;
    }
  }
  return (e / 2.0);
}

//----------------------------------------------------------------------------//

void rpmd_pile::init(size_t ndim, size_t natom, size_t nbead, const double &kT,
                     const double &dt, const double *mass,
                     const double *cartpos, const double *cartvel,
                     double gamma_centroid) {
  // Intialize the RPMD base class
  rpmd_base::init(ndim, natom, nbead, kT, dt, mass, cartpos, cartvel,
                  gamma_centroid);

  c1 = arma::vec(nbead);
  c2 = arma::vec(nbead);

  init_langevin(dt, gamma_centroid);
}

//----------------------------------------------------------------------------//

void rpmd_pile::init_langevin(const double &dt, double gamma_centroid) {
  for (size_t k = 0; k < nbeads(); ++k) {
    double gamma(0);
    if (k == 0) {
      gamma = gamma_centroid;
    } else {
      gamma = 2.0 * m_omega_k(k);
    }

    c1(k) = std::exp(-0.5 * dt * gamma);
    c2(k) = std::sqrt(1.0 - c1(k) * c1(k));
  }

  std::random_device rd{};
  pile_prng.seed(rd());
}

void rpmd_pile::seed_pile_prng(int my_seed) { pile_prng.seed(my_seed); }

//----------------------------------------------------------------------------//

void rpmd_pile::step(const double &dt) {
  const double sqrt_beta_n = std::sqrt(m_beta_n);

  // 1. Advance thermostats dt2
  mom_c2n();
  for (size_t i = 0; i < ndofs(); ++i) {
    double fac = m_sqrt_mass(i) / sqrt_beta_n;
    for (size_t k = 0; k < nbeads(); ++k) {
      double rand_gaussian = rand_01(pile_prng);
      m_mom_nmode(i, k) =
          m_mom_nmode(i, k) * c1(k) + fac * c2(k) * rand_gaussian;
    }
  }
  mom_n2c();

  // 2. Do un-thermostatted RPMD evolution
  rpmd_base::step(dt);

  // 3. Advance thermostats final dt2 and recalc KE
  mom_c2n();
  for (size_t i = 0; i < ndofs(); ++i) {
    double fac = m_sqrt_mass(i) / sqrt_beta_n;
    for (size_t k = 0; k < nbeads(); ++k) {
      double rand_gaussian = rand_01(pile_prng);
      m_mom_nmode(i, k) =
          m_mom_nmode(i, k) * c1(k) + fac * c2(k) * rand_gaussian;
    }
  }

  mom_n2c();
  m_Ekin = calc_KE();

  m_temp_kT = m_Ekin * 2.0 / ndofs() / nbeads(); // not actual temperature, kT
}

//----------------------------------------------------------------------------//

double rpmd_pile::invariant() const { return m_Ekin + m_Espring + m_Epot_sum; }

//----------------------------------------------------------------------------//

} // namespace parts

////////////////////////////////////////////////////////////////////////////////
