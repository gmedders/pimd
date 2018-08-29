#include <random>

#include "eom.h"

namespace parts {

////////////////////////////////////////////////////////////////////////////////

void eom::print_params() {
  m_potential->print_params();
  std::cout << "# gammaTh_fac = " << m_gamma << std::endl;
}

//----------------------------------------------------------------------------//

void eom::set_gammaTh(const double &dt, double gam_fac) {
  m_gamma = gam_fac * m_potential->get_w();
  // FIXME
  // init_langevin(dt, m_gamma);
}

//----------------------------------------------------------------------------//

double eom::force(size_t ndim, size_t natom, size_t nbeads, const double *x,
                  double *f) {
  return m_potential->force(ndim, natom, nbeads, x, f);
}

//----------------------------------------------------------------------------//

void eom::set_up(const size_t ndim, const size_t natom, const size_t nbead,
                 const double beta, const double dt, double *pos, double *vel) {
  size_t ndof = ndim * natom;

  // Generate initial positions, if needed
  std::vector<double> all_bead_crd;
  if (pos == nullptr) {
    for (size_t i = 0; i < nbead * ndof; ++i) {
      all_bead_crd.push_back(0.0);
    }
    pos = &all_bead_crd[0];
  }

  // Generate initial velocities, if needed
  std::vector<double> all_bead_vel;
  if (vel == nullptr) {
    double kT = 1.0 / beta;

    std::random_device rd{};
    std::mt19937 prng(rd());
    std::normal_distribution<> rand_01{0, 1};

    for (size_t i = 0; i < nbead * ndof; ++i) {
      const double sigma = std::sqrt(kT / atm_mass);
      double v = sigma * rand_01(prng);
      all_bead_vel.push_back(v);
    }
    vel = &all_bead_vel[0];
  }

  // Now populate the mass array
  double mass[ndof];
  for (size_t i = 0; i < natom; ++i) {
    const double Mi = atm_mass;
    for (size_t k = 0; k < ndim; ++k)
      mass[k + ndim * i] = Mi;
  }

  // setup the simulation
  m_potential->set_params(params);

  // FIXME
  // // PIMD
  // init(ndim, natom, nbead, 1.0 / beta, mass, pos, vel);
  // // RPMD
  // init(ndim, natom, nbead, 1.0 / beta, dt, mass, pos, vel, m_gamma);
  // // VV
  // init(ndim, natom, dt, mass, pos, vel);
}

////////////////////////////////////////////////////////////////////////////////

} // namespace parts

////////////////////////////////////////////////////////////////////////////////
