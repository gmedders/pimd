#include <cassert>
#include <cmath>

#include <algorithm>

#include "vv_base.h"

////////////////////////////////////////////////////////////////////////////////

namespace parts {

//----------------------------------------------------------------------------//

void vv_base::init(size_t ndim, size_t natom, const double &dt,
                   const double *mass, const double *cartpos,
                   const double *cartvel) {
  assert(ndim * natom > 0);
  assert(mass != 0 && cartpos != 0 && cartvel != 0);

  m_ndim = ndim;
  m_natoms = natom;
  m_ndof = ndim * natom;
  m_dt = dt;

  m_pos = arma::vec(m_ndof);
  m_mom = arma::vec(m_ndof);
  m_frc = arma::vec(m_ndof);

  m_mass = arma::vec(m_ndof);

  // initialize cartesian positions and momenta
  for (size_t i = 0; i < m_ndof; ++i) {
    m_pos(i) = cartpos[i];
    m_mom(i) = mass[i] * cartvel[i];
  }

  for (size_t i = 0; i < m_ndof; ++i) {
    m_mass(i) = mass[i];
  }

  // calculate the initial kinetic energy

  m_Ekin = 0.0;
  for (size_t i = 0; i < m_ndof; ++i) {
    m_Ekin += m_mom(i) * m_mom(i) / m_mass(i);
  }

  m_Ekin /= 2;

  m_Epot = force(m_ndim, m_natoms, 1, m_pos.memptr(), m_frc.memptr());
}

//----------------------------------------------------------------------------//

void vv_base::step(const double &dt) {
  const double dt2 = dt / 2;

  // Following equations 21-25 of dx.doi.org/10.1063/1.3489925
  // 1. Evolution of RP momenta under Hamiltonian V_{n}[q(t0)] by dt/2
  m_mom += dt2 * m_frc;

  // 2. Advance positions
  m_pos += dt * m_mom / m_mass;

  // 3. Final evolution of RP momenta under Hamiltonian V_{n}[q(t0+dt)]
  m_Epot = force(m_ndim, m_natoms, 1, m_pos.memptr(), m_frc.memptr());

  m_mom += dt2 * m_frc;

  m_Ekin = 0.0;
  for (size_t i = 0; i < m_ndof; ++i) {
    m_Ekin += m_mom(i) * m_mom(i) / m_mass(i);
  }

  m_Ekin /= 2;
  m_temp_kT = m_Ekin * 2.0 / m_ndof; // not actual temperature, kT
}

//----------------------------------------------------------------------------//

double vv_base::invariant() const {
  // Epot is in kcal/mol already
  return m_Ekin + m_Epot;
}

//----------------------------------------------------------------------------//

} // namespace parts

////////////////////////////////////////////////////////////////////////////////
