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

#include "helpers.h"

#include "sim-classes.h"

#define DUMP_TRAJ yes

//
// constant temperature PIMD
// units are au
//

////////////////////////////////////////////////////////////////////////////////

namespace {

//const size_t print_time = 10000; // au
//const size_t equil_time = 10000000;

//test
const size_t print_time = 1; // au
const size_t equil_time = 0;

} // namespace

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
    // 1. load the coordinates

    std::cout.setf(std::ios_base::showpoint);
    std::cout.precision(9);

    if (argc != 5) {
        std::cerr << "usage: classical-md beta dt gammaTh_factor nsamples" << std::endl;
        return EXIT_FAILURE;
    }

    size_t ndim = 2;
    size_t natom = 2;
    size_t nbead = 1;

    double beta = parts::parse_to_double(argv[1]);
    double dt = parts::parse_to_double(argv[2]);
    double gammaTh_fac = parts::parse_to_double(argv[3]);
    int nsamples = parts::parse_to_int(argv[4]);

    if(nsamples <= 0)
        nsamples = 1000;

    const size_t prod_time = print_time*nsamples;
    const size_t simulation_time = equil_time + prod_time; // au

    const size_t nsteps = int(simulation_time / dt);
    const size_t nsteps_equil = int(equil_time / dt);
    const size_t nprint = int(print_time / dt);

    //rpmd sim;
    parts::rpmd sim;
    sim.m_potential->set_all_bead_states(0, nbead);

    //sim.m_potential->set_bath_params(ndim, natom - 1,
    //                               sim.m_gamma, 2*0.9*sim.m_potential->get_w(),
    //                                sim.m_potential->get_m());

    try {
        int nx=4;
        double x[] = {1.75, -2.5, -1.75, -4.5};
        double v[] = {0.001, 0.001, -0.001, 0.001};

        std::vector<double> all_crd;
        std::vector<double> all_vel;
        int count(0);
        for(int i = 0; i < nbead*ndim*natom; ++i){
            //all_vel.push_back(0.0);
            if(count == nx)
                count = 0;

            all_crd.push_back(x[count]);
            all_vel.push_back(v[count]);
            ++count;
        }

        sim.set_up(ndim, natom, nbead, beta, dt, &all_crd[0], &all_vel[0]);
    } catch (const std::exception& e) {
        std::cerr << " ** Error ** : " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    sim.set_gammaTh(dt, gammaTh_fac);
    sim.print_params();

    // 2. iterate
#ifdef DUMP_TRAJ
    std::ostringstream ss_filename;
    ss_filename << "cart_traj-class_" << int(beta) << ".dat";

    std::ofstream of_cart_traj;
    of_cart_traj.open(ss_filename.str().c_str());
    of_cart_traj << std::scientific;
    of_cart_traj.precision(15);
#endif

    int count(0);
    double sum_Espring(0);
    for (size_t n = 0; n < nsteps; ++n) {
        sim.step(dt);
        //if (n%nprint == 0) {
        if (n%nprint == 0 && n >= nsteps_equil) {
            ++count;
            sum_Espring += sim.Espring();
            std::cout << n*dt << ' '
                      << sim.invariant() << ' '
                      << sim.m_potential->avg_active_state() << ' '
                      << sim.Ep() << ' '
                      << sim.temp_kT() << ' '
                      << sim.avg_cart_pos() << std::endl;

#ifdef DUMP_TRAJ
            sim.dump_1D_frame(of_cart_traj);
#endif
        }

    }

#ifdef DUMP_TRAJ
    of_cart_traj.close();
#endif

    return EXIT_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
