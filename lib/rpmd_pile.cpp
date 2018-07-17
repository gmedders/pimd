#include <cassert>
#include <cmath>

#include <algorithm>

#include "rpmd_pile.h"

#define DECOMPOSE_KE yes

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

void rpmd_pile::init(size_t ndof, size_t nbead, const double &kT,
                     const double &dt, const double *mass,
                     const double *cartpos, const double *cartvel,
                     double gamma_centroid) {
  // Intialize the RPMD base class
  rpmd_base::init(ndof, nbead, kT, dt, mass, cartpos, cartvel, gamma_centroid);

  // Now do all the RPMD-NHC specific stuff

  c1 = arma::vec(nbead);
  c2 = arma::vec(nbead);

  // rand_gaussian = arma::mat(ndof, nbead, arma::zeros);

  for (size_t k = 0; k < nbead; ++k) {
    double gamma(0);
    if (k == 0) {
      gamma = gamma_centroid;
    } else {
      gamma = 2.0 * m_omega_k(k);
    }

    c1(k) = std::exp(-0.5 * dt * gamma);
    c2(k) = std::sqrt(1.0 - c1(k) * c1(k));
  }

#ifdef DECOMPOSE_KE
  saved_mom = arma::mat(ndof, nbead);
#endif

  std::random_device rd{};
  pile_prng.seed(rd());

  // std::cerr << "<<< Thermostatting ( tau = " << 1.0/gamma_centroid
  //          << " ) >>>"<<std::endl;
}

void rpmd_pile::seed_pile_prng(int my_seed) { pile_prng.seed(my_seed); }

//----------------------------------------------------------------------------//

void rpmd_pile::step(const double &dt) {
  const double dt2 = dt / 2;
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

#ifdef DECOMPOSE_KE
  // Store the correct nmode momenta
  saved_mom = m_mom_nmode;

  // Zero non-centroid normal mode momenta
  for (size_t n = 1; n < nbeads(); ++n) {
    for (size_t i = 0; i < ndofs(); ++i) {
      m_mom_nmode(i, n) = 0.0;
    }
  }
  mom_n2c();
  double Ekin_centroid = calc_KE();

  // Now restore m_mom_nmode and zero-out centoid momenta
  m_mom_nmode = saved_mom;
  for (size_t n = 0; n < 1; ++n) {
    for (size_t i = 0; i < ndofs(); ++i) {
      m_mom_nmode(i, n) = 0.0;
    }
  }
  mom_n2c();
  double Ekin_higherNM = calc_KE();

  // Finally, restore m_mom_nmode and calculate full KE
  m_mom_nmode = saved_mom;
#endif
  mom_n2c();
  m_Ekin = calc_KE();
#ifndef DECOMPOSE_KE
  double Ekin_centroid = m_Ekin;
  double Ekin_higherNM = m_Ekin;
#endif

  m_temp_kT = m_Ekin * 2.0 / ndofs() / nbeads(); // not actual temperature, kT
  m_temp_kT_centroid = Ekin_centroid * 2.0 / ndofs();
  if (nbeads() > 1) {
    m_temp_kT_higherNM = Ekin_higherNM * 2.0 / ndofs() / (nbeads() - 1);
  } else {
    m_temp_kT_higherNM = m_temp_kT;
  }
}

//----------------------------------------------------------------------------//

double rpmd_pile::invariant() const { return m_Ekin + m_Espring + m_Epot_sum; }

//----------------------------------------------------------------------------//

} // namespace parts

////////////////////////////////////////////////////////////////////////////////
