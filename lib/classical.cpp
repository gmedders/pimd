#include <cassert>
#include <cmath>

#include <algorithm>

#include "classical.h"

////////////////////////////////////////////////////////////////////////////////

namespace parts {

//----------------------------------------------------------------------------//

void classical::init(size_t ndim, size_t natom, const double &dt,
                     const double beta, const double atm_mass,
                     const double *cartpos, const double *cartvel) {
  // Set m_ndim, m_natom, m_ndof, m_dt
  eom::set_up(ndim, natom, dt);

  assert(ndim * natom > 0);

  m_pos = arma::vec(m_ndof);
  m_vel = arma::vec(m_ndof);
  m_mom = arma::vec(m_ndof);
  m_frc = arma::vec(m_ndof);

  m_mass = arma::vec(m_ndof);

  if (cartpos == nullptr) {
    m_pos.zeros();
  } else {
    std::copy(cartpos, cartpos + m_ndof, m_pos.memptr());
  }

  if (cartvel == nullptr) {
    generate_boltzmann_velocities(beta, atm_mass, m_vel.memptr());
  } else {
    std::copy(cartvel, cartvel + m_ndof, m_vel.memptr());
  }

  std::fill(m_mass.memptr(), m_mass.memptr() + m_ndof, atm_mass);

  // Having initialized velocities and masses, set the momenta
  m_mom = m_vel * m_mass;

  // calculate the initial kinetic energy
  m_Ekin = 0.0;
  for (size_t i = 0; i < m_ndof; ++i) {
    m_Ekin += m_mom(i) * m_mom(i) / m_mass(i);
  }
  m_Ekin /= 2;

  m_Epot = force(m_ndim, m_natom, 1, m_pos.memptr(), m_frc.memptr());
} // namespace parts

void classical::calc_pos_stats(void) {
  assert(ndim() * natoms() == 1);

  m_avg_cart_pos = m_pos[0];
}

void classical::dump_1D_frame(std::ostream &of_traj) {
  of_traj << m_potential->state_id[0] << ' ';
  for (size_t i = 0; i < ndofs(); ++i) {
    of_traj << ' ' << m_pos(i) << ' ' << m_mom(i) / m_mass(i);
  }
  of_traj << std::endl;
}

//----------------------------------------------------------------------------//

} // namespace parts

////////////////////////////////////////////////////////////////////////////////
