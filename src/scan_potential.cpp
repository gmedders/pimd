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

#include "sim-classes.h"

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
  // 1. load the coordinates

  std::cout.setf(std::ios_base::showpoint);
  std::cout.precision(9);

  double dt = 1.0;
  double beta = 20.0;

  int nbead(1);
  int ndim(2);
  int natom(1);

  // rpmd sim;
  parts::vv sim;
  sim.m_potential->set_all_bead_states(1, nbead);
  double hop_params[] = {0.02, dt, beta, 0.0};
  sim.m_potential->set_hopping_params(hop_params);

  double crd[] = {0.0, 0.0};
  double frc[] = {0.0, 0.0};

  try {
    sim.set_up(ndim, natom, nbead, beta, dt, crd, frc);
  } catch (const std::exception &e) {
    std::cerr << " ** Error ** : " << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  for (double x = -100.0; x < 100.0; x += 0.1) {
    sim.cart_ptr()[0] = x;
    std::cout << x << ' ';
    sim.m_potential->set_all_bead_states(0, nbead);
    std::cout << sim.force(ndim, natom, nbead, sim.cart_ptr(), frc) << ' ';
    sim.m_potential->set_all_bead_states(1, nbead);
    std::cout << sim.force(ndim, natom, nbead, sim.cart_ptr(), frc) << ' ';
    std::cout << std::endl;
    // sim.m_pos(0) = x;
    // sim.force(1, 1, sim.m_pos.memptr(), frc);
  }

  return EXIT_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
