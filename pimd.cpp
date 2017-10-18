#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif // HAVE_CONFIG_H

#include <cmath>
#include <cassert>
#include <cstdlib>

#include <sstream>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "nhc.h"
#include "helpers.h"

#include "sim-classes.h"

//
// constant temperature PIMD
// units are au
//

////////////////////////////////////////////////////////////////////////////////

namespace {

//const size_t print_time = 100; // au
//const size_t equil_time = 100000;
//const size_t prod_time = 100000;
const size_t print_time = 100; // au
const size_t equil_time = 5000000;
const size_t prod_time = 10000000;
const size_t simulation_time = equil_time + prod_time; // au

} // namespace

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
    // 1. load the coordinates

    std::cout.setf(std::ios_base::showpoint);
    std::cout.precision(9);

    if (argc != 4) {
        std::cerr << "usage: pimd nbeads beta dt" << std::endl;
        return EXIT_FAILURE;
    }

    size_t ndim = 1;
    size_t natom = 1;

    size_t nbead = strtod(argv[1], NULL);
    double beta = parts::parse_to_double(argv[2]);
    double dt = parts::parse_to_double(argv[3]);

    const size_t nsteps = int(simulation_time / dt);
    const size_t nsteps_equil = int(equil_time / dt);
    const size_t nprint = int(print_time / dt);

    //pimd sim;
    parts::pimd sim;
    sim.m_potential.set_all_bead_states(0, nbead);
    double GammaEl(0.0);
    double voltage(0.0);
    double hop_params[] = {GammaEl, dt, beta, voltage};
    sim.m_potential.set_hopping_params(hop_params);

    try {
        sim.set_up(nbead, ndim, natom, beta);
    } catch (const std::exception& e) {
        std::cerr << " ** Error ** : " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    // 2. iterate
    std::ostringstream ss_filename;
    ss_filename << "cart_traj-pimd_" << nbead << '_' << int(beta) << ".dat";

    std::ofstream of_cart_traj;
    of_cart_traj.open(ss_filename.str().c_str());
    of_cart_traj << std::scientific;
    of_cart_traj.precision(15);

    for (size_t n = 0; n < nsteps; ++n) {
        sim.step(dt);
        if (n%nprint == 0 && n >= nsteps_equil) {
            std::cout << n*dt << ' '
                      << sim.invariant() << ' '
                      << sim.Espring() << ' '
                      << sim.Ek() << ' '
                      << sim.Ep() << ' '
                      << sim.temp_kT() << ' '
                      << sim.avg_cart_pos() << std::endl;

            sim.dump_1D_frame(of_cart_traj);
        }

    }

    of_cart_traj.close();

    return EXIT_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
