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

typedef pot::sho potential_type;

const size_t nsteps = 20000; 

double atm_mass(2000); // au

const double dt = 1; // au

const size_t nframe = 1000; // 0.5 ps
const size_t nprint = 10;

////////////////////////////////////////////////////////////////////////////////

struct pimd : public parts::pimd_base {

    ~pimd();

    void set_up(const size_t, const size_t, const size_t, const double);
    double force(const double*, double*);

    inline size_t natom() const { return m_natom; }
    inline double Ep() const { return m_Epot_sum; }

    void set_COM_at_origin();

private:
    size_t m_natom;
    size_t m_ndim;
    size_t m_ndofs;
    double* pos;
    potential_type m_potential;
};

//----------------------------------------------------------------------------//

pimd::~pimd()
{
    delete[] pos;
}

void pimd::set_up(const size_t nbead, const size_t ndim, const size_t natom,
                  const double beta)
{
    m_natom = natom;
    m_ndim = ndim;

    m_ndofs = m_natom*m_ndim;

    pos = new double[m_ndofs];

    // prepare masses
    double* mass = new double[m_ndofs]; // for every degree of freedom
    for (size_t i = 0; i < natom; ++i) {
        const double Mi = atm_mass;
        for (size_t k = 0; k < ndim; ++k)
            mass[k + ndim*i] = Mi;
    }

    // setup the simulation

    init(m_ndofs, nbead, 1.0/beta, mass, pos);

    // clean up

    delete[] mass;
}

//----------------------------------------------------------------------------//

double pimd::force(const double* x, double* f)
{
    double Epot = m_potential(x[0], f[0]);

//    for (size_t n = 0; n < ndim*natom(); ++n)
//        f[n] *= -engunit;

    return Epot;
}

////////////////////////////////////////////////////////////////////////////////

} // namespace

int main(int argc, char** argv)
{
    // 1. load the coordinates
    
    std::cout.setf(std::ios_base::showpoint);
    std::cout.precision(9);

    if (argc != 3) {
        std::cerr << "usage: pimd nbeads T" << std::endl;
        return EXIT_FAILURE;
    }

    size_t nbead = 8; 
    size_t ndim = 3;
    size_t natom = 1;

    nbead = strtod(argv[1], NULL);

    double beta(1.0);
    {
        std::istringstream iss(argv[3]);
        iss >> beta;
        if (!iss || !iss.eof()) {
            std::cerr << "could not convert '" << argv[2]
                      << "' to real number" << std::endl;
            return EXIT_FAILURE;
        }

        assert(beta > 0.0);
    }

    pimd sim;

    try {
        sim.set_up(nbead, ndim, natom, beta);
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

    }

    return EXIT_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
