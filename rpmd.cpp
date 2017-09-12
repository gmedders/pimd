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
#include "mt19937.h"

#include "sim_classes.h"

//
// constant temperature PIMD
// units are au
//

////////////////////////////////////////////////////////////////////////////////

namespace {

const size_t print_time = 1; // au
const size_t equil_time = 0;
const size_t prod_time = 1000;

//const size_t print_time = 20; // au
//const size_t equil_time = 1000;
//const size_t prod_time = 100000;
const size_t simulation_time = equil_time + prod_time; // au

} // namespace

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
    // 1. load the coordinates

    std::cout.setf(std::ios_base::showpoint);
    std::cout.precision(9);

    if (argc != 4) {
        std::cerr << "usage: rpmd nbeads beta dt" << std::endl;
        return EXIT_FAILURE;
    }

    size_t ndim = 1;
    size_t natom = 1;

    size_t nbead = strtod(argv[1], NULL);

    double beta(1.0);
    {
        std::istringstream iss(argv[2]);
        iss >> beta;
        if (!iss || !iss.eof()) {
            std::cerr << "could not convert '" << argv[2]
                      << "' to real number" << std::endl;
            return EXIT_FAILURE;
        }

        assert(beta > 0.0);
    }

    double dt(1.0);
    {
        std::istringstream iss(argv[3]);
        iss >> dt;
        if (!iss || !iss.eof()) {
            std::cerr << "could not convert '" << argv[3]
                      << "' to real number" << std::endl;
            return EXIT_FAILURE;
        }

        assert(dt > 0.0);
    }

    const size_t nsteps = int(simulation_time / dt);
    const size_t nsteps_equil = int(equil_time / dt);
    const size_t nprint = int(print_time / dt);

    //rpmd sim;
    parts::rpmd sim;

    try {
        sim.set_up_new_init_cond(nbead, ndim, natom, beta, dt);
    } catch (const std::exception& e) {
        std::cerr << " ** Error ** : " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    // 2. iterate
    std::ofstream of_cart_traj;
    of_cart_traj.open("cart_traj.dat");
    of_cart_traj << std::scientific;
    of_cart_traj.precision(15);

    for (size_t n = 0; n < nsteps; ++n) {
        sim.step(dt);
        if (n%nprint == 0) {
        //if (n%nprint == 0 && n >= nsteps_equil) {
            std::cout << n*dt << ' '
                      << sim.invariant() << ' '
                      << sim.Espring() << ' '
                      << sim.Ek() << ' '
                      << sim.Ep() << ' '
                      << sim.temp_kT() << ' '
                      << sim.avg_cart_pos() << std::endl;

            //sim.dump_1D_frame(of_cart_traj);
        }

    }

    of_cart_traj.close();

    return EXIT_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
