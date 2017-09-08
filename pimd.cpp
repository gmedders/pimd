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

#include "pimd-base.h"

#include "sho.h"

//
// constant temperature PIMD
// units are au
//

namespace {

////////////////////////////////////////////////////////////////////////////////

typedef pot::sho potential_type;

double atm_mass(2000); // au
double omega(0.2); // omega

const size_t print_time = 20; // au
const size_t equil_time = 100000;
const size_t prod_time = 100000;
const size_t simulation_time = equil_time + prod_time; // au

////////////////////////////////////////////////////////////////////////////////

struct pimd : public parts::pimd_base {

    ~pimd();

    void set_up(const size_t, const size_t, const size_t, const double);
    double force(const double*, double*);

    inline double Espring() const { return m_Espring; }
    inline double Ep() const { return m_Epot_sum; }
    inline double Ek() const { return m_Ekin_fict; }
    inline double temp_kT() const { return m_temp_kT; }
    double avg_cart_pos(void);

    void dump_1D_frame(std::ofstream&);

private:
    size_t m_natom;
    size_t m_ndim;
    size_t m_ndofs;
    size_t m_nbead;
    double* pos;
    potential_type m_potential;
};

//----------------------------------------------------------------------------//

pimd::~pimd()
{
    delete[] pos;
}

//----------------------------------------------------------------------------//

double pimd::avg_cart_pos(void)
{
    double avg(0);
    for(size_t i = 0; i < m_nbead; ++i)
        avg += m_pos_cart[i*m_ndofs];
    return avg/m_nbead;
}

//----------------------------------------------------------------------------//

void pimd::dump_1D_frame(std::ofstream& of_traj)
{
    vel_n2c();

    double ekin_cart(0);
    of_traj << m_nbead << std::endl;
    for(size_t i = 0; i < m_nbead; ++i) {
        for(size_t j = 0; j < m_ndofs; ++j) {
            of_traj << ' ' << m_pos_cart[i*m_ndofs + j]
                    << ' ' << m_vel_cart[i*m_ndofs + j];
            ekin_cart += atm_mass
                         * m_vel_cart[i*m_ndofs + j]*m_vel_cart[i*m_ndofs + j];
        }
        of_traj << std::endl;
    }
#if 0
//    of_traj << "m_pos_cart:  ";
//    for(size_t i = 0; i < m_nbead; ++i){
//        for(size_t j = 0; j < m_ndofs; ++j){
//            of_traj << ' ' << m_pos_cart[i*m_ndofs + j];
//        }
//    }
//    of_traj << std::endl;
//
//    of_traj << "m_vel_nmode: ";
//    for(size_t i = 0; i < m_nbead; ++i){
//        for(size_t j = 0; j < m_ndofs; ++j){
//            of_traj << ' ' << m_vel_nmode[i*m_ndofs + j];
//        }
//    }
//    of_traj << std::endl;
//
//    of_traj << "m_vel_cart:  ";
//    for(size_t i = 0; i < m_nbead; ++i){
//        for(size_t j = 0; j < m_ndofs; ++j){
//            of_traj << ' ' << m_vel_cart[i*m_ndofs + j];
//            m_vel_nmode[i*m_ndofs + j] = 0.0;
//            ekin_cart += atm_mass
//                         * m_vel_cart[i*m_ndofs + j]*m_vel_cart[i*m_ndofs + j];
//        }
//    }
//    of_traj << std::endl;
//
//    vel_c2n();
//
//    of_traj << "new_v_nmode: ";
//    for(size_t i = 0; i < m_nbead; ++i){
//        for(size_t j = 0; j < m_ndofs; ++j){
//            of_traj << ' ' << m_vel_nmode[i*m_ndofs + j];
//        }
//    }
//    of_traj << std::endl;
#endif

    ekin_cart /= 2.0;
//    of_traj << " Ekin_nm= " << m_Ekin_fict
//            << " Ekin_cart= " << ekin_cart << std::endl;
}

//----------------------------------------------------------------------------//

void pimd::set_up(const size_t nbead, const size_t ndim, const size_t natom,
                  const double beta)
{
    m_natom = natom;
    m_ndim = ndim;
    m_nbead = nbead;

    m_ndofs = m_natom*m_ndim;

    pos = new double[m_ndofs];

    pos[0] = 0.0;

    // prepare masses
    double* mass = new double[m_ndofs]; // for every degree of freedom
    for (size_t i = 0; i < natom; ++i) {
        const double Mi = atm_mass;
        for (size_t k = 0; k < ndim; ++k)
            mass[k + ndim*i] = Mi;
    }

    // setup the simulation

    m_potential.set_params(omega, atm_mass);

    init(m_ndofs, nbead, 1.0/beta, mass, pos);

    // clean up

    delete[] mass;
}

//----------------------------------------------------------------------------//

double pimd::force(const double* x, double* f)
{
    double Epot = m_potential(m_natom, x, f);

    return Epot;
}

////////////////////////////////////////////////////////////////////////////////

} // namespace

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

    pimd sim;

    try {
        sim.set_up(nbead, ndim, natom, beta);
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

            sim.dump_1D_frame(of_cart_traj);
        }

    }

    of_cart_traj.close();

    return EXIT_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
