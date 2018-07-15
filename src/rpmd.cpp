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

#include "helpers.h"

#include "sim_classes.h"

//
// constant temperature PIMD
// units are au
//

////////////////////////////////////////////////////////////////////////////////

namespace {

const size_t print_time = 10000; // au
const size_t equil_time = 10000000;

// test
// const size_t print_time = 1; // au
// const size_t equil_time = 0;

} // namespace

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
  // 1. load the coordinates

  std::cout.setf(std::ios_base::showpoint);
  std::cout.precision(9);

  if (argc != 8) {
    std::cerr
        << "usage: rpmd nbeads beta dt GammaEl gammaTh_factor voltage nsamples"
        << std::endl;
    return EXIT_FAILURE;
  }

  size_t ndim = 1;
  size_t natom = 1;

  int nbead = parts::parse_to_int(argv[1]);
  double beta = parts::parse_to_double(argv[2]);
  double dt = parts::parse_to_double(argv[3]);
  double GammaEl = parts::parse_to_double(argv[4]);
  double gammaTh_fac = parts::parse_to_double(argv[5]);
  double voltage = parts::parse_to_double(argv[6]);
  int nsamples = parts::parse_to_int(argv[7]);

  if (nsamples <= 0)
    nsamples = 1000;

  const size_t prod_time = print_time * nsamples;
  const size_t simulation_time = equil_time + prod_time; // au

  const size_t nsteps = int(simulation_time / dt);
  const size_t nsteps_equil = int(equil_time / dt);
  const size_t nprint = int(print_time / dt);

  // rpmd sim;
  parts::rpmd sim;
  sim.m_potential.set_all_bead_states(0, nbead);
  double hop_params[] = {GammaEl, dt, beta / nbead, voltage};
  sim.m_potential.set_hopping_params(hop_params);
  sim.set_gammaTh(gammaTh_fac);

  try {
    // sim.set_up_new_init_cond(nbead, ndim, natom, beta, dt);
    // 64 bead example starting configuration. take slices from it to seed
    // lower-number beads
    int nx = 64;
    double x[] = {
        1.229025298970351e+01, 1.228969606844854e+01, 1.228565692653722e+01,
        1.228395801534676e+01, 1.229303568829728e+01, 1.229020485881905e+01,
        1.228835139115360e+01, 1.228838234168862e+01, 1.229142269596170e+01,
        1.229378215208285e+01, 1.229110520094307e+01, 1.229112515078028e+01,
        1.229793576108822e+01, 1.229848335179882e+01, 1.229307236160982e+01,
        1.229200974081691e+01, 1.229025298970351e+01, 1.228969606844854e+01,
        1.228565692653722e+01, 1.228395801534676e+01, 1.229303568829728e+01,
        1.229020485881905e+01, 1.228835139115360e+01, 1.228838234168862e+01,
        1.229142269596170e+01, 1.229378215208285e+01, 1.229110520094307e+01,
        1.229112515078028e+01, 1.229793576108822e+01, 1.229848335179882e+01,
        1.229307236160982e+01, 1.229200974081691e+01, 1.229025298970351e+01,
        1.228969606844854e+01, 1.228565692653722e+01, 1.228395801534676e+01,
        1.229303568829728e+01, 1.229020485881905e+01, 1.228835139115360e+01,
        1.228838234168862e+01, 1.229142269596170e+01, 1.229378215208285e+01,
        1.229110520094307e+01, 1.229112515078028e+01, 1.229793576108822e+01,
        1.229848335179882e+01, 1.229307236160982e+01, 1.229200974081691e+01,
        1.229025298970351e+01, 1.228969606844854e+01, 1.228565692653722e+01,
        1.228395801534676e+01, 1.229303568829728e+01, 1.229020485881905e+01,
        1.228835139115360e+01, 1.228838234168862e+01, 1.229142269596170e+01,
        1.229378215208285e+01, 1.229110520094307e+01, 1.229112515078028e+01,
        1.229793576108822e+01, 1.229848335179882e+01, 1.229307236160982e+01,
        1.229200974081691e+01};
    std::vector<double> all_crd;
    std::vector<double> all_vel;
    int count(0);
    for (int i = 0; i < nbead * ndim * natom; ++i) {
      all_vel.push_back(0.0);
      if (count == nx)
        count = 0;

      all_crd.push_back(x[count]);
      ++count;
    }
    // Center the ring polymer

    double avg(0);
    for (int i = 0; i < nbead; ++i) {
      avg += all_crd[i];
    }
    avg /= nbead;

    for (int i = 0; i < nbead; ++i) {
      all_crd[i] -= avg + 1.6;
    }

    sim.set_up(nbead, ndim, natom, beta, dt, &all_crd[0], &all_vel[0]);
  } catch (const std::exception &e) {
    std::cerr << " ** Error ** : " << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  sim.print_params();

  // 2. iterate
  std::ostringstream ss_filename;
  ss_filename << "cart_traj-rpmd_" << nbead << '_' << int(beta) << ".dat";

  std::ofstream of_cart_traj;
  of_cart_traj.open(ss_filename.str().c_str());
  of_cart_traj << std::scientific;
  of_cart_traj.precision(15);

  int count(0);
  double sum_Espring(0);
  for (size_t n = 0; n < nsteps; ++n) {
    sim.step(dt);
    // if (n%nprint == 0) {
    if (n % nprint == 0 && n >= nsteps_equil) {
      ++count;
      sum_Espring += sim.Espring();
      std::cout << n * dt << ' ' << sim.invariant() << ' ' << sim.Espring()
                << ' ' << sum_Espring / count
                << ' '
                //<< sim.Ek() << ' '
                << sim.m_potential.avg_active_state() << ' ' << sim.Ep() << ' '
                << sim.temp_kT() << ' ' << sim.temp_kT_centroid() << ' '
                << sim.temp_kT_higherNM() << ' ' << sim.avg_cart_pos()
                << std::endl;

      sim.dump_1D_frame(of_cart_traj);
    }
  }

  of_cart_traj.close();

  return EXIT_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
