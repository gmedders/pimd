#include <cassert>
#include <cmath>

#include <algorithm>
#include <random>

#include "rpmd_nhc.h"

#include "nhc.h"

////////////////////////////////////////////////////////////////////////////////

namespace parts {

//----------------------------------------------------------------------------//

rpmd_nhc::rpmd_nhc() : rpmd_base() { m_thermostats = 0; }

//----------------------------------------------------------------------------//

rpmd_nhc::~rpmd_nhc() { delete[] m_thermostats; }

//----------------------------------------------------------------------------//

void rpmd_nhc::init(size_t ndim, size_t natom, size_t nbead, const double &kT,
                    const double &dt, const double *mass, const double *cartpos,
                    const double *cartvel, double dummy) {

  // Intialize the RPMD base class
  rpmd_base::init(ndim, natom, nbead, kT, dt, mass, cartpos, cartvel, dummy);
  size_t ndof = ndim * natom;

  // Now do all the RPMD-NHC specific stuff

  std::cerr << "<<< Thermostatting >>>" << std::endl;
  delete[] m_thermostats;

  const size_t th_size = nhc::size(nchain);
  m_thermostats = new double[ndof * nbead * th_size];

  std::mt19937 prng;
  prng.seed(19107);

  m_omega_M = kT / hbar;
  m_tau = 2 * M_PI / m_omega_M;

  for (size_t b = 0; b < nbead; ++b)
    for (size_t i = 0; i < ndof; ++i) {
      const size_t j = i + b * ndof;
      nhc::initialize(nchain, m_thermostats + th_size * j, m_tau, prng);
    }
}

//----------------------------------------------------------------------------//

void rpmd_nhc::step(const double &dt) {
  const double dt2 = dt / 2;

  const size_t th_size = nhc::size(nchain);

  // 1. Advance thermostats dt2 and rescale velocities
  for (size_t n = 0; n < nbeads(); ++n) {
    for (size_t i = 0; i < ndofs(); ++i) {
      const size_t j = n * ndofs() + i;
      const double mass = m_mass(i);
      const double Ekin2 = m_mom_cart(i, n) * m_mom_cart(i, n) / mass;

      const double aa = nhc::advance(nchain, m_thermostats + j * th_size, m_tau,
                                     Ekin2 / m_kT, dt2);
      m_mom_cart(i, n) *= aa;
    }
  }

  // 2. Do un-thermostatted RPMD evolution
  rpmd_base::step(dt);

  // 3. Advance thermostats final dt2, rescale velocities, and recalc KE
  m_Ekin = 0.0;
  for (size_t n = 0; n < nbeads(); ++n) {
    for (size_t i = 0; i < ndofs(); ++i) {
      double mass = m_mass(i);
      const double Ekin2 = m_mom_cart(i, n) * m_mom_cart(i, n) / mass;
      const size_t j = n * ndofs() + i;
      const double aa = nhc::advance(nchain, m_thermostats + j * th_size, m_tau,
                                     Ekin2 / m_kT, dt2);

      m_mom_cart(i, n) *= aa;
      m_Ekin += Ekin2 * aa * aa;
    }
  }

  m_Ekin /= 2;
  m_temp_kT = m_Ekin * 2.0 / ndofs() / nbeads(); // not actual temperature, kT
}

//----------------------------------------------------------------------------//

double rpmd_nhc::invariant() const {
  const size_t th_size = nhc::size(nchain);

  double accu(0);
  for (size_t n = 0; n < nbeads() * ndofs(); ++n)
    accu += nhc::invariant(nchain, m_thermostats + n * th_size, m_tau);

  // Epot is in kcal/mol already
  return (m_Ekin + m_Espring + m_kT * accu) / engunit + m_Epot_sum;
}

//----------------------------------------------------------------------------//

} // namespace parts

////////////////////////////////////////////////////////////////////////////////
