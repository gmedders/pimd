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

const double print_time = 0.05/0.0002; // au
const double prod_time = 60.0/0.0002; // au
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

    std::cout << "# w  = " << parts::omega << std::endl;
    std::cout << "# m  = " << parts::atm_mass << std::endl;
    std::cout << "# g  = " << parts::bb_x0 << std::endl;
    std::cout << "# dG = " << parts::dG << std::endl;

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

    std::vector<double> sum_pos;
    std::vector<double> sum_pos_L1;
    std::vector<double> sum_pos_L2;
    std::vector<double> sum_pos_Linf;
    std::vector<double> sum_temp;
    std::vector<double> sum_state;

    std::vector<double> traj_pos;
    std::vector<double> traj_pos_L1;
    std::vector<double> traj_pos_L2;
    std::vector<double> traj_pos_Linf;
    std::vector<double> traj_temp;
    std::vector<double> traj_sum_state;


    for (size_t i = 0; i < nsteps; ++i) {
        double itime = i*dt;
        if (i%nprint == 0 && itime <= tcf_max_time) {
            nsamples.push_back(0);
            tcf.push_back(0.0);
            time.push_back(itime);

            sum_pos.push_back(0.0);
            sum_pos_L1.push_back(0.0);
            sum_pos_L2.push_back(0.0);
            sum_pos_Linf.push_back(0.0);
            sum_temp.push_back(0.0);
            sum_state.push_back(0.0);

            traj_pos.push_back(0.0);
            traj_pos_L1.push_back(0.0);
            traj_pos_L2.push_back(0.0);
            traj_pos_Linf.push_back(0.0);
            traj_temp.push_back(0.0);
            traj_sum_state.push_back(0.0);
        }
    }
    const size_t tcf_max_nsteps = time.size();

    int ntemp(0);

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
        std::vector<int> init_active_state;

        for(size_t n = 0; n < nbead; ++n){
            std::string line;
            std::getline(ifs, line);
            ++lineno;
            std::istringstream iss(line);

            int this_state;
            iss >> this_state;
            init_active_state.push_back(this_state);

            for(size_t i = 0; i < ndof; ++i){
                double q;
                double v;
                iss >> q >> v;
                all_bead_crd.push_back(q);
                all_bead_vel.push_back(v);
            }
            check_parsing(iss, lineno);
        }

        //std::cerr << ntemp << ' ' << std::endl;

        // Now set up this simulation
        
        //rpmd sim;
        //parts::vv sim;
        parts::rpmd sim;
        sim.m_potential.set_individual_bead_states(init_active_state);
        double GammaEl(1.0e-3);
        double hop_params[] = {GammaEl, dt, beta};
        sim.m_potential.set_hopping_params(hop_params);

        //std::ostringstream ss_filename;
        //ss_filename << "traj_" << ntemp << ".dat";
        //std::ofstream of_cart_traj;
        //of_cart_traj.open(ss_filename.str());
        //of_cart_traj << std::scientific;
        //of_cart_traj.precision(15);

        try {
            //sim.set_up_new_init_cond(nbead, ndim, natom, beta, dt,
            //                         &all_bead_crd[0]);
            sim.set_up(nbead, ndim, natom, beta, dt,
                       &all_bead_crd[0], &all_bead_vel[0]);
        } catch (const std::exception& e) {
            std::cerr << " ** Error ** : " << e.what() << std::endl;
            return EXIT_FAILURE;
        }

        size_t count(0);
        for (size_t n = 0; n < nsteps; ++n) {
            sim.step(dt);
            if (n%nprint == 0) {
                sim.calc_pos_stats();
                traj_pos[count] = sim.avg_cart_pos();
                traj_pos_L1[count] = sim.L1_cart_pos();
                traj_pos_L2[count] = sim.L2_cart_pos();
                traj_pos_Linf[count] = sim.Linf_cart_pos();
                traj_temp[count] = sim.temp_kT(); // kT
                traj_sum_state[count] = sim.m_potential.sum_active_state();
                ++count;
            }
        }

#if 0
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
#endif
        // Accumulate the Average Temperature
        for (size_t i = 0; i < sum_temp.size(); ++i){
            sum_pos[i] += traj_pos[i];
            sum_pos_L1[i] += traj_pos_L1[i];
            sum_pos_L2[i] += traj_pos_L2[i];
            sum_pos_Linf[i] += traj_pos_Linf[i];
            sum_temp[i] += traj_temp[i];
            sum_state[i] += traj_sum_state[i]/nbead;
        }
        ++ntemp;
    }

    // Finally, print the results
    std::cout << std::scientific;
    std::cout.precision(10);
    for (size_t i = 0; i < tcf.size(); ++i){
        std::cout << std::setw(20) << time[i]
//                  << std::setw(20) << tcf[i]/nsamples[i]
                  << std::setw(20) << sum_state[i]/ntemp
                  << std::setw(20) << sum_temp[i]/ntemp
                  << std::setw(20) << sum_pos[i]/ntemp
                  << std::setw(20) << sum_pos_L1[i]/ntemp
                  << std::setw(20) << sum_pos_L2[i]/ntemp
                  << std::setw(20) << sum_pos_Linf[i]/ntemp
                  << std::endl;
    }


    ifs.close();

    return EXIT_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
