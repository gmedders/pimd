#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif // HAVE_CONFIG_H

#ifdef ENABLE_MPI
#include <mpi.h>
#endif

#include <cmath>
#include <cassert>
#include <cstdlib>

#include <sstream>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "helpers.h"

#include "sim-classes.h"

////////////////////////////////////////////////////////////////////////////////

namespace {

int my_rank(0), my_size(1);

//const double print_time = 1.0/0.0002; // au
//const double prod_time = 10000.0/0.0002; // au

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

#ifdef ENABLE_MPI
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &my_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif

    //int rand_seed[] = {1507040009, 1507040010, 1507040011, 1507040012};
    //srand(rand_seed[my_rank]);
    srand(time(NULL) + my_rank);

    if (argc < 4) {
        if(my_rank == 0)
            std::cerr << "usage: ensemble_tcf_rpmd input_file dt GammaEl time"
                      << std::endl;
        return EXIT_FAILURE;
    }

    size_t ndim = 1;
    size_t natom = 1;

    double dt = parts::parse_to_double(argv[2]);
    double GammaEl = parts::parse_to_double(argv[3]);


    double prod_time;

    if(argc == 5){
        prod_time = parts::parse_to_double(argv[4]);
    } else {
        prod_time = 200.0/0.0002; // au
    }

    const double print_time = prod_time/5000; // au
    const double tcf_max_time = prod_time;

    size_t nsteps = int(prod_time / dt);
    size_t nprint = int(print_time / dt);

    if(my_rank == 0){
        std::cout << "# w  = " << parts::omega << std::endl;
        std::cout << "# m  = " << parts::atm_mass << std::endl;
        std::cout << "# g  = " << parts::bb_x0 << std::endl;
        std::cout << "# dG = " << parts::dG << std::endl;
    }

    // 2. iterate
    std::string filename(argv[1]);
    std::ifstream ifs(filename.c_str());

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
    std::vector<double> traj_temp_count;
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
            traj_temp_count.push_back(0.0);
            traj_sum_state.push_back(0.0);
        }
    }
    const size_t tcf_max_nsteps = time.size();

    int iframe(0);

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

        // Every process reads the entire file, but only run the simulation
        // if the frame id matches process rank
        if(iframe % my_size != my_rank){
            ++iframe;
            continue;
        }else{
            ++iframe;
        }


#ifdef ENABLE_MPI
        //std::cerr << "Rank " << my_rank << " running frame " << iframe << std::endl;
#endif

        // Now set up this simulation
        
        //rpmd sim;
        //parts::vv sim;
        parts::rpmd sim;
        sim.m_potential.set_individual_bead_states(init_active_state);
        double hop_params[] = {GammaEl, dt, beta};
        sim.m_potential.set_hopping_params(hop_params);

        try {
            //sim.set_up_new_init_cond(nbead, ndim, natom, beta, dt,
            //                         &all_bead_crd[0]);
            sim.set_up(nbead, ndim, natom, beta, dt,
                       &all_bead_crd[0], &all_bead_vel[0]);
        } catch (const std::exception& e) {
            std::cerr << " ** Error ** : " << e.what() << std::endl;
            return EXIT_FAILURE;
        }

        //std::fill(traj_temp_count.begin(), traj_temp_count.end(), 0.0);

        //std::ostringstream ss_filename;
        //ss_filename << "traj_" << iframe << "_proc" << my_rank << ".dat";
        //std::ofstream of_cart_traj;
        //of_cart_traj.open(ss_filename.str());
        //of_cart_traj << std::scientific;
        //of_cart_traj.precision(15);

        //of_cart_traj << "# " << all_bead_crd[0] << ' ' << all_bead_vel[0]
        //             << std::endl;

        size_t count(0);
        for (size_t n = 0; n < nsteps; ++n) {
            sim.step(dt);
            if (n%nprint == 0) {
                //sim.calc_pos_stats();
                traj_pos[count] = sim.avg_cart_pos();
                traj_pos_L1[count] = sim.L1_cart_pos();
                traj_pos_L2[count] = sim.L2_cart_pos();
                traj_pos_Linf[count] = sim.Linf_cart_pos();
                traj_sum_state[count] = sim.m_potential.sum_active_state();

                //of_cart_traj << time[count] << ' ' << traj_sum_state[count]
                //    << ' ' << sim.m_potential.prev_rand() <<std::endl;

                //if(sim.m_potential.sum_active_state() == 0){
                    traj_temp[count] = sim.temp_kT(); // kT
                    traj_temp_count[count] += 1.0;
                //}
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
            //std::cerr << time[i] << ' ' << sum_state[i]/iframe << std::endl;
        }
    }

#ifdef ENABLE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &sum_state[0], sum_state.size(),
                  MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &sum_temp[0], sum_temp.size(),
                  MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &sum_pos[0], sum_pos.size(),
                  MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &sum_pos_L1[0], sum_pos_L1.size(),
                  MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &sum_pos_L2[0], sum_pos_L2.size(),
                  MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &sum_pos_Linf[0], sum_pos_Linf.size(),
                  MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &traj_temp_count[0], traj_temp_count.size(),
                  MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
#endif

    if(my_rank == 0){
        // Finally, print the results
        std::cout << std::scientific;
        std::cout.precision(10);
        for (size_t i = 0; i < tcf.size(); ++i){
            std::cout << std::setw(20) << time[i]
                //                  << std::setw(20) << tcf[i]/nsamples[i]
                << std::setw(20) << sum_state[i]/iframe
                << std::setw(20) << sum_temp[i]/traj_temp_count[i]
                << std::setw(20) << sum_pos[i]/iframe
                << std::setw(20) << sum_pos_L1[i]/iframe
                << std::setw(20) << sum_pos_L2[i]/iframe
                << std::setw(20) << sum_pos_Linf[i]/iframe
                << std::endl;
        }
    }

    ifs.close();

#ifdef ENABLE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
#endif

    return EXIT_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
