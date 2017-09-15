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

#include "sim-classes.h"

////////////////////////////////////////////////////////////////////////////////

namespace {

const double print_time = 1.0; // au
const double prod_time = 60/0.003; // au
//const double prod_time = 200;

const double tcf_max_time = prod_time;
//const double tcf_max_time = 30;
const double simulation_time = prod_time; // au

void check_parsing(std::istringstream& iss, size_t lineno)
{
    if (iss.fail()) {
        std::ostringstream oss;
        oss << "failed to parse line " << lineno
            << " of the input file:" << std::endl << iss.str() << std::endl;
        throw std::runtime_error(oss.str());
    }
}

} // namespace

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
    // 1. load the coordinates

    std::cout.setf(std::ios_base::showpoint);
    std::cout.precision(9);

    double dt = 1.0;
    double beta = 20.0;

    int nbead(1);
    int ndim(1);
    int natom(1);


    //rpmd sim;
    parts::vv sim;
    sim.m_potential.set_active_state(1);
    double hop_params[] = {0.02, dt, beta};
    sim.m_potential.set_hopping_params(hop_params);

    double crd[1] = {0.0};
    double frc[1] = {0.0};

    try {
        sim.set_up(nbead, ndim, natom, beta, dt,
                   crd, frc);
    } catch (const std::exception& e) {
        std::cerr << " ** Error ** : " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    for(double x = -100.0; x < 100.0; x += 0.1){
        sim.cart_ptr()[0] = x;
        std::cout << x << ' ';
        sim.force(1, 1, sim.cart_ptr(), frc);
        //sim.m_pos(0) = x;
        //sim.force(1, 1, sim.m_pos.memptr(), frc);
    }

    return EXIT_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
