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

#include "sim_classes.h"

////////////////////////////////////////////////////////////////////////////////

namespace {

const double print_time = 0.1; // au
const double prod_time = 200;

const double tcf_max_time = 30;
const double simulation_time = prod_time; // au

void check_parsing(std::istringstream& iss, size_t lineno)
{
    if (iss.fail()) {
        std::ostringstream oss;
        oss << "failed to parse line " << lineno
            << " of the input file:" << std::endl << iss << std::endl;
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

    if (argc != 3) {
        std::cerr << "usage: ensemble_tcf_rpmd input_file dt" << std::endl;
        return EXIT_FAILURE;
    }

    size_t ndim = 1;
    size_t natom = 1;

    double dt(1.0);
    {
        std::istringstream iss(argv[2]);
        iss >> dt;
        if (!iss || !iss.eof()) {
            std::cerr << "could not convert '" << argv[2]
                      << "' to real number" << std::endl;
            return EXIT_FAILURE;
        }

        assert(dt > 0.0);
    }

    const size_t nsteps = int(simulation_time / dt);
    const size_t nprint = int(print_time / dt);

    // 2. iterate
    std::string filename(argv[1]);
    std::ifstream ifs(filename);

    size_t lineno(0);

    // Initialize arrays
    std::vector<double> nsamples;
    std::vector<double> tcf;
    std::vector<double> time;

    for (size_t i = 0; i < nsteps; ++i) {
        double itime = i*dt;
        if (i%nprint == 0 && itime <= tcf_max_time) {
            nsamples.push_back(0);
            tcf.push_back(0.0);
            time.push_back(i*dt);
        }
    }
    const size_t tcf_max_nsteps = time.size();

    while(!ifs.eof()){

        // Read this frame of the input file
        // Fist line: NBead NDof Beta
        std::string line;
        std::getline(ifs, line);
        ++lineno;

        if(ifs.eof())
            break;

        size_t nbead;
        size_t ndof;
        double beta;
        std::istringstream iss(line);
        iss >> nbead >> ndof >> beta;
        check_parsing(iss, lineno);

        assert(nbead > 0);
        assert(ndof == ndim*natom);
        assert(beta > 0);

        // Next NBead lines: q1 v1 q2 v2
        std::vector<double> all_bead_crd;
        std::vector<double> all_bead_vel;

        for(size_t n = 0; n < nbead; ++n){
            std::string line;
            std::getline(ifs, line);
            ++lineno;
            std::istringstream iss(line);

            for(size_t i = 0; i < ndof; ++i){
                double q;
                double v;
                iss >> q >> v;
                check_parsing(iss, lineno);
                all_bead_crd.push_back(q);
                all_bead_vel.push_back(v);
            }
        }

        // Now set up this simulation
        
        //rpmd sim;
        parts::rpmd sim;

        try {
            sim.set_up(nbead, ndim, natom, beta, dt,
                       &all_bead_crd[0], &all_bead_vel[0]);
        } catch (const std::exception& e) {
            std::cerr << " ** Error ** : " << e.what() << std::endl;
            return EXIT_FAILURE;
        }

        std::vector<double> traj_pos;

        for (size_t n = 0; n < nsteps; ++n) {
            sim.step(dt);
            if (n%nprint == 0) {
                //std::cout << n*dt << ' '
                //          << sim.invariant() << ' '
                //          << sim.Espring() << ' '
                //          << sim.Ek() << ' '
                //          << sim.Ep() << ' '
                //          << sim.temp_kT() << ' '
                //          << sim.avg_cart_pos() << std::endl;
                traj_pos.push_back(sim.avg_cart_pos());
            }
        }

        // Accumulate the TCF
        for (size_t i = 0; i < traj_pos.size(); ++i){
            for (size_t j = i; (j < traj_pos.size())
                            && ((j - i) < tcf_max_nsteps); ++j)
            {

                size_t delta = j - i;
                nsamples[delta] += 1.0;

                tcf[delta] += traj_pos[i] * traj_pos[j];
            }
        }
    }

    // Finally, print the results
    std::cout << std::scientific;
    std::cout.precision(10);
    for (size_t i = 0; i < tcf.size(); ++i){
        std::cout << std::setw(20) << time[i]
                  << std::setw(20) << tcf[i]/nsamples[i]
                  << std::endl;
    }


    ifs.close();

    return EXIT_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
