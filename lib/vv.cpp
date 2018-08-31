#include <cassert>
#include <cmath>

#include <algorithm>

#include "vv.h"

////////////////////////////////////////////////////////////////////////////////

namespace parts {

//----------------------------------------------------------------------------//

void vv::step(const double &dt) {
  const double dt2 = dt / 2;

  // Following equations 21-25 of dx.doi.org/10.1063/1.3489925
  // 1. Evolution of RP momenta under Hamiltonian V_{n}[q(t0)] by dt/2
  m_mom += dt2 * m_frc;

  // 2. Advance positions
  m_pos += dt * m_mom / m_mass;

  // 3. Final evolution of RP momenta under Hamiltonian V_{n}[q(t0+dt)]
  m_Epot = force(m_ndim, m_natom, 1, m_pos.memptr(), m_frc.memptr());

  m_mom += dt2 * m_frc;

  m_Ekin = 0.0;
  for (size_t i = 0; i < m_ndof; ++i) {
    m_Ekin += m_mom(i) * m_mom(i) / m_mass(i);
  }

  m_Ekin /= 2;
  m_temp_kT = m_Ekin * 2.0 / m_ndof; // not actual temperature, kT
}

//----------------------------------------------------------------------------//

double vv::invariant() const {
  // Epot is in kcal/mol already
  return m_Ekin + m_Epot;
}

//----------------------------------------------------------------------------//

} // namespace parts

////////////////////////////////////////////////////////////////////////////////
