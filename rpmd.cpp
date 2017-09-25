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

#include "sim-classes.h"

//
// constant temperature PIMD
// units are au
//

////////////////////////////////////////////////////////////////////////////////

namespace {

const size_t print_time = 10000; // au
//10000
//const size_t equil_time = 100000000;
//const size_t prod_time = 100000000;
//1000
const size_t equil_time = 10000000;
const size_t prod_time = 10000000;
//test
//const size_t equil_time = 0;
////const size_t prod_time = 30000000;
//const size_t prod_time = 20*0.0002;

//const size_t print_time = 1000; // au
//const size_t equil_time = 0;
//const size_t prod_time = 5000;

//const size_t print_time = 1; // au
//const size_t equil_time = 0;
//const size_t prod_time = 10000000;
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
    sim.m_potential.set_active_state(0);
    double GammaEl(0.0);
    double hop_params[] = {GammaEl, dt, beta};
    sim.m_potential.set_hopping_params(hop_params);

    try {
        //sim.set_up_new_init_cond(nbead, ndim, natom, beta, dt);
        // 32 bead example starting configuration. take slices from it to seed lower-number beads
        //double x[] = { 0.0};
        double x[] = { 1.229025298970351e+01,1.228969606844854e+01,1.228565692653722e+01,1.228395801534676e+01,1.229303568829728e+01,1.229020485881905e+01,1.228835139115360e+01,1.228838234168862e+01,1.229142269596170e+01,1.229378215208285e+01,1.229110520094307e+01,1.229112515078028e+01,1.229793576108822e+01,1.229848335179882e+01,1.229307236160982e+01,1.229200974081691e+01,1.229025298970351e+01,1.228969606844854e+01,1.228565692653722e+01,1.228395801534676e+01,1.229303568829728e+01,1.229020485881905e+01,1.228835139115360e+01,1.228838234168862e+01,1.229142269596170e+01,1.229378215208285e+01,1.229110520094307e+01,1.229112515078028e+01,1.229793576108822e+01,1.229848335179882e+01,1.229307236160982e+01,1.229200974081691e+01};
        double v[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        sim.set_up(nbead, ndim, natom, beta, dt, x, v);
    } catch (const std::exception& e) {
        std::cerr << " ** Error ** : " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    // 2. iterate
    std::ostringstream ss_filename;
    ss_filename << "cart_traj-rpmd_" << nbead << '_' << int(beta) << ".dat";

    std::ofstream of_cart_traj;
    of_cart_traj.open(ss_filename.str());
    of_cart_traj << std::scientific;
    of_cart_traj.precision(15);

    for (size_t n = 0; n < nsteps; ++n) {
        sim.step(dt);
        //if (n%nprint == 0) {
        if (n%nprint == 0 && n >= nsteps_equil) {
            std::cout << n*dt << ' '
                      << sim.invariant() << ' '
                      << sim.Espring() << ' '
                      << sim.Ek()*beta << ' '
                      //<< sim.Ek() << ' '
                      << sim.m_potential.active_state << ' '
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
