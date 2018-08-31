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

void eom::set_up(const size_t ndim, const size_t natom, const double dt) {
  m_ndim = ndim;
  m_natom = natom;
  m_ndof = ndim * natom;

  m_dt = dt;

  // // Now populate the mass array
  // double mass[ndof];
  // for (size_t i = 0; i < natom; ++i) {
  //   const double Mi = atm_mass;
  //   for (size_t k = 0; k < ndim; ++k)
  //     mass[k + ndim * i] = Mi;
  // }

  m_potential->set_params(params);
}

void eom::generate_boltzmann_velocities(const double beta,
                                        const double atm_mass, double *vel) {

  double kT = 1.0 / beta;

  std::random_device rd{};
  std::mt19937 prng(rd());
  std::normal_distribution<> rand_01{0, 1};

  for (size_t i = 0; i < m_ndof; ++i) {
    const double sigma = std::sqrt(kT / atm_mass);
    double v = sigma * rand_01(prng);
    vel[i] = v;
  }
}

////////////////////////////////////////////////////////////////////////////////

} // namespace parts

////////////////////////////////////////////////////////////////////////////////
