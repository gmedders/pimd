#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <cassert>
#include <cmath>
#include <cstdlib>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "gcmc.h"

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
  // 1. load the coordinates

  std::cout.setf(std::ios_base::showpoint);
  std::cout.precision(4);

  const size_t ndim(2);
  const size_t natom(2);
  const size_t nbead(1);

  double beta(1024);
  double dt(1.0e-12);

  double all_bead_crd[] = {1.0, 1.1, 2.0, 2.1};
  double all_bead_vel[] = {0.5, 0.6, 0.7, 0.8};

  int active_state(0);

  const size_t nsteps = 2;
  const size_t nprint = 1;

  // 2. Set up GCMC simulation object
  parts::gcmc sim;

  sim.set_up(ndim, natom, nbead, beta, dt, all_bead_crd, all_bead_vel);

  sim.m_md_ensemble->m_potential->set_all_bead_states(0, nbead);
  // double GammaEl(1.0e-3);
  // double hop_params[] = {GammaEl, dt, beta, 0.0};
  // sim.m_md_ensemble->m_potential->set_hopping_params(hop_params);

  sim.dump(std::cout);

  sim.step(dt, beta);

  sim.dump(std::cout);

  return EXIT_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
