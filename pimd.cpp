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

#include "read-xyz.h"

#include "pimd-base.h"

//
// constant temperature PIMD
// units are au
//

namespace {

typedef x2o::pot_water_ion potential_type;

const size_t nsteps = 20000; 
size_t nbead = 8; 
size_t ndim = 3;

double atm_mass(2000); // au

static double  T_kelvin = 200.0; // [Kelvin]
const double dt = 1; // au

const size_t nframe = 1000; // 0.5 ps
const size_t nprint = 10;

////////////////////////////////////////////////////////////////////////////////

struct pimd : public parts::pimd_base {

    void load(const char* filename);
    double force(const double*, double*);

    inline size_t natom() const { return m_natom; }
    inline double Ep() const { return m_Epot_sum; }

    void set_COM_at_origin();
    void append_frame(const char* filename, const double& time);

private:
    size_t m_natom;
    potential_type m_potential;
};

//----------------------------------------------------------------------------//

void pimd::load(const char* filename)
{

    // Load the configuration

    //std::vector<tools::xyz_frame> molecs;

    //try {
    //    tools::read_xyz(filename, molecs);
    //} catch (const std::exception& e) {
    //    std::cerr << " ** Error ** : " << e.what() << std::endl;
    //    exit(1);
    //}

    double* pos = molecs[0].xyz;

    const size_t natom = 1;
    m_natom = natom;

    // prepare masses
    double* mass = new double[ndim*natom]; // for every degree of freedom
    for (size_t i = 0; i < natom; ++i) {
        const double Mi = atm_mass;
        for (size_t k = 0; k < ndim; ++k)
            mass[k + ndim*i] = Mi;
    }

    // setup the simulation

    init(ndim*natom, nbead, kB*T, mass, pos);

    // clean up

    delete[] mass;
}

//----------------------------------------------------------------------------//

void pimd::set_COM_at_origin()
{
    double com[ndim] = {0.0}, mass_sum(0);

    // centroid's COM
    for (size_t n = 0; n < ndim*natom(); n += ndim) {
        mass_sum += m_fict_mass[n];
        for (size_t k = 0; k < ndim; ++k)
            com[k] += m_fict_mass[n]*m_pos_nmode[n + k];
    }

    for (size_t k = 0; k < ndim; ++k)
        com[k] /= mass_sum;

    for (size_t n = 0; n < ndim*natom(); n += ndim)
        for (size_t k = 0; k < ndim; ++k)
            m_pos_nmode[n + k] -= com[k];

    pos_n2c();
}

//----------------------------------------------------------------------------//

double pimd::force(const double* x, double* f)
{
    double Epot = m_potential(x, f);

    for (size_t n = 0; n < ndim*natom(); ++n)
        f[n] *= -engunit;

    return Epot;
}

////////////////////////////////////////////////////////////////////////////////

} // namespace

int main(int argc, char** argv)
{
    // 1. load the coordinates
    
    std::cout.setf(std::ios_base::showpoint);
    std::cout.precision(9);

    if (argc != 5) {
        std::cerr << "usage: pimd nbeads in.xyz T out.xyz" << std::endl;
        return EXIT_FAILURE;
    }

    nbead = strtod(argv[1], NULL);

    const char* input_filename = argv[2];
    const char* output_filename = argv[4];

    {
        std::istringstream iss(argv[3]);
        iss >> T;
        if (!iss || !iss.eof()) {
            std::cerr << "could not convert '" << argv[2]
                      << "' to real number" << std::endl;
            return EXIT_FAILURE;
        }

        assert(T > 0.0);
    }

    pimd sim;

    try {
        sim.load(input_filename);
    } catch (const std::exception& e) {
        std::cerr << " ** Error ** : " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    std::cerr << "molecule loaded!" << std::endl
	      << "  > natom = " << sim.natom() << std::endl ;

    // 2. iterate

    for (size_t n = 0; n < nsteps; ++n) {
        sim.step(dt);
        if (n%nprint == 0) {
            std::cout << n*dt << ' '
                      << sim.invariant() << ' '
                      << sim.Ep() << std::endl;
        }

        if (n%nframe == 0)
            sim.append_frame(output_filename, n*dt);
    }

    return EXIT_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////

void pimd::append_frame(const char* filename, const double& t)
{
    set_COM_at_origin();

    std::ofstream ofs(filename, std::ios_base::app);
    ofs << natom()*nbeads() << '\n'
        << "t = " << t << ", nbeads = " << nbeads()
        << ", T = " << T << '\n';

    for (size_t b = 0; b < nbeads(); ++b){
        for (size_t w = 0; w < natom(); ++w) {
            print_xyz(ofs, 'X', m_pos_cart + b*ndofs() + ndim*i);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
