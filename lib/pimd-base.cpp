#include <cassert>
#include <cmath>

#include "mt19937.h"

#include "pimd-base.h"

#define DO_NHC yes

#ifdef DO_NHC
#include "nhc.h"
#endif

////////////////////////////////////////////////////////////////////////////////

namespace parts {

//----------------------------------------------------------------------------//

pimd_base::pimd_base() : necklace() {
  m_fict_mass = 0;
  m_thermostats = 0;
}

//----------------------------------------------------------------------------//

pimd_base::~pimd_base() {
  delete[] m_fict_mass;
#ifdef DO_NHC
  delete[] m_thermostats;
#endif
}

//----------------------------------------------------------------------------//

void pimd_base::init(size_t ndim, size_t natom, size_t nbead, const double &kT,
                     const double *mass, const double *cartpos) {
  assert(ndim * natom > 0 && nbead > 0);
  assert(nbead % 2 == 0 || nbead == 1);
  assert(kT > 0.0 && mass != 0 && cartpos != 0);

  size_t ndof = ndim * natom;

  necklace::setup(ndim, natom, nbead);

  delete[] m_fict_mass;
#ifdef DO_NHC
  delete[] m_thermostats;
#endif

  m_kT = kT;
  m_omega_M = std::sqrt(double(nbead)) * kT / hbar;

  // fictitious masses

  m_fict_mass = new double[ndof * nbead];

  for (size_t b = 0; b < nbead; ++b) {
    const double factor = (b == 0 ? 1.0 : lambda(b));
    for (size_t i = 0; i < ndof; ++i)
      m_fict_mass[i + b * ndof] = mass[i] * factor;
  }

  // thermostats

#ifdef DO_NHC
  const size_t th_size = nhc::size(nchain);
  m_thermostats = new double[ndof * nbead * th_size];
#endif

  mt19937 prg(27606);

  m_tau = 2 * M_PI / m_omega_M;

#ifdef DO_NHC
  for (size_t b = 0; b < nbead; ++b)
    for (size_t i = 0; i < ndof; ++i) {
      const size_t j = i + b * ndof;
      nhc::initialize(nchain, m_thermostats + th_size * j, m_tau, prg);
    }
#endif

  // initialize cartesian positions

  for (size_t b = 0; b < nbead; ++b)
    for (size_t i = 0; i < ndof; ++i)
      m_pos_cart[i + b * ndof] = cartpos[i];

  pos_c2n(); // cart -> nmode

  // initialize normal mode velocities

  m_Ekin_fict = 0.0;
  for (size_t b = 0; b < nbead; ++b)
    for (size_t i = 0; i < ndof; ++i) {
      const size_t j = i + b * ndof;
      const double sigma = std::sqrt(kT / m_fict_mass[j]);
      m_vel_nmode[j] = sigma * prg.random_gaussian();
      m_Ekin_fict += m_fict_mass[j] * m_vel_nmode[j] * m_vel_nmode[j];
    }

  m_Ekin_fict /= 2;

  //// scale kinetic energy to the target temperature

  // const double scale_factor = std::sqrt(0.5*kT*nbead*ndof/m_Ekin_fict);
  // for (size_t b = 0; b < nbead; ++b)
  //    for (size_t i = 0; i < ndof; ++i)
  //        m_vel_nmode[i + b*ndof] *= scale_factor;

  // m_Ekin_fict *= scale_factor*scale_factor;

  // compute forces

  pimd_force();
}

//----------------------------------------------------------------------------//

void pimd_base::pimd_force() {
  const size_t n_total = ndofs() * nbeads();

  // zero out m_frc_cart
  std::fill(m_frc_cart, m_frc_cart + n_total, 0.0);

  // compute forces for each bead
  m_Epot_sum = force(ndim(), natoms(), nbeads(), m_pos_cart, m_frc_cart);

  // For PIMD, divide energy by number of beads
  m_Epot_sum /= nbeads();

  // transform it to normal modes
  frc_c2n();

  const double omega2 = m_omega_M * m_omega_M;

  m_Espring = 0.0;

  // add the harmonic part

  for (size_t b = 1; b < nbeads(); ++b) {
    for (size_t i = 0; i < ndofs(); ++i) {
      const size_t j = b * ndofs() + i;
      const double tmp = m_fict_mass[j] * omega2 * m_pos_nmode[j];
      m_frc_nmode[j] -= tmp;
      m_Espring += tmp * m_pos_nmode[j];
    }
  }

  m_Espring /= 2;
}

//----------------------------------------------------------------------------//

void pimd_base::step(const double &dt) {
#ifdef DO_NHC
  const size_t th_size = nhc::size(nchain);
#endif

  // 1. advance thermostats, velocities by dt/2, nmode position on dt

  const double dt2 = dt / 2;

  for (size_t b = 0; b < nbeads(); ++b) {
    for (size_t i = 0; i < ndofs(); ++i) {
      const size_t j = b * ndofs() + i;
      const double mass = m_fict_mass[j];
      const double Ekin2 = mass * m_vel_nmode[j] * m_vel_nmode[j];

#ifdef DO_NHC
      const double aa = nhc::advance(nchain, m_thermostats + j * th_size, m_tau,
                                     Ekin2 / m_kT, dt2);
#else
      const double aa = 1.0;
#endif
      m_vel_nmode[j] = aa * m_vel_nmode[j] + dt2 * m_frc_nmode[j] / mass;
      m_pos_nmode[j] += dt * m_vel_nmode[j];
    }
  }

  // 2. transform position to cartesian, compute forces

  pos_n2c();

  pimd_force(); // computes normal mode forces

  // 3. advance velocities and thermostats by dt/2

  m_Ekin_fict = 0.0;

  for (size_t b = 0; b < nbeads(); ++b) {
    for (size_t i = 0; i < ndofs(); ++i) {
      const size_t j = b * ndofs() + i;
      const double mass = m_fict_mass[j];

      m_vel_nmode[j] += dt2 * m_frc_nmode[j] / mass;
      const double Ekin2 = mass * m_vel_nmode[j] * m_vel_nmode[j];
#ifdef DO_NHC
      const double aa = nhc::advance(nchain, m_thermostats + j * th_size, m_tau,
                                     Ekin2 / m_kT, dt2);

      m_vel_nmode[j] *= aa;
      m_Ekin_fict += Ekin2 * aa * aa;
#else
      m_Ekin_fict += Ekin2;
#endif
    }
  }

  m_Ekin_fict /= 2;
  m_temp_kT =
      m_Ekin_fict * 2.0 / ndofs() / nbeads(); // not actual temperature, kT
}

//----------------------------------------------------------------------------//

double pimd_base::invariant() const {
  double accu(0);

#ifdef DO_NHC
  const size_t th_size = nhc::size(nchain);

  for (size_t n = 0; n < nbeads() * ndofs(); ++n)
    accu += nhc::invariant(nchain, m_thermostats + n * th_size, m_tau);
#endif

  // Epot is in kcal/mol already
  return (m_Ekin_fict + m_Espring + m_kT * accu) / engunit + m_Epot_sum;
}

//----------------------------------------------------------------------------//

} // namespace parts

////////////////////////////////////////////////////////////////////////////////
