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
#include "constraining-sphere.h"

//#include "mbpol.h"
#include "pot-water-ion.h"

//
// constant temperature PIMD (with "constraining sphere")
// units are like in DLPOLY: length = A,
// time = ps, and energy conversion factor
// is in pimd_base.h
//

namespace {

typedef x2o::pot_water_ion potential_type;
//typedef ttm::ttm3f potential_type;
//typedef ttm::ps::pot_nasa potential_type;

const size_t nsteps = 20000; 
size_t nbead = 8; 

static double  T = 200.0; // [degrees]
const double dt = 0.0001; // time-step [ps]

const size_t nframe = 1000; // 0.5 ps
const size_t nprint = 10;

const double H_mass = 1.0079;
const double O_mass = 15.9949;
const double I_mass = 126;

double X_params[] = { -1 ,
                     65.60240487*0.14818471,
		     4155.38,
		     3.32398,
		     31778.2,
		     2.70971,
		     892.645,
		     1.78065,
		     4944.06,
		     2.74377};

////////////////////////////////////////////////////////////////////////////////

void print_xyz(std::ostream& os, char atom, const double* x)
{
    os << std::setw(2) << std::left << atom
       << std::setprecision(9) << std::scientific
       << std::setw(18) << std::right << x[0]
       << std::setw(18) << std::right << x[1]
       << std::setw(18) << std::right << x[2]
       << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

struct pimd : public parts::pimd_base {

    void load(const char* filename);
    double force(const double*, double*);

    inline size_t nw()    const { return m_nw; }
    inline size_t nI()    const { return m_nI; }
    inline size_t natom() const { return m_natom; }
    inline double Ep() const { return m_Epot_sum; }

    double Es() const;

    void set_COM_at_origin();
    void append_frame(const char* filename, const double& time);

private:
    size_t m_nw;
    size_t m_nI;
    size_t m_natom;
    potential_type m_potential;
    parts::constraining_sphere m_sphere;
};

//----------------------------------------------------------------------------//

void pimd::load(const char* filename)
{

    // Load the configuration

    std::vector<tools::xyz_frame> molecs;

    try {
	tools::read_xyz(filename, molecs);
    } catch (const std::exception& e) {
	std::cerr << " ** Error ** : " << e.what() << std::endl;
	exit(1);
    }

    size_t nw = molecs[0].natm/3;
    size_t nI = 1;
    double* pos = molecs[0].xyz;

    const size_t natom = 3*nw + nI;
    m_nw = nw;
    m_nI = nI;
    m_natom = natom;

    // define potential
    m_potential.define_ion_params(X_params);

    // prepare masses

    double* mass = new double[3*natom]; // for every degree of freedom
    for (size_t i = 0; i < 3*m_nw; ++i) {
        const double Mi = (i%3 == 0 ? O_mass : H_mass);
        for (size_t k = 0; k < 3; ++k)
            mass[k + 3*i] = Mi;
    }
    for (size_t i = 0; i < nI; ++i) {
        const double Mi = I_mass;
     	for (size_t k = 0; k < 3; ++k)
            mass[k + 3*(3*m_nw + i)] = Mi;
    }

    // setup the simulation

    init(3*natom, nbead, kB*T, mass, pos);

    m_sphere.setup(4.5, 100.0);

    // clean up

    delete[] mass;
}

//----------------------------------------------------------------------------//

double pimd::Es() const
{
    double E(0);

    for (size_t b = 0; b < nbeads(); ++b)
        E += m_sphere(m_nw, m_pos_cart + b*ndofs(), 0);

    return E;
}

//----------------------------------------------------------------------------//

void pimd::set_COM_at_origin()
{
    double com[3] = {0.0, 0.0, 0.0}, mass_sum(0);

    // centroid's COM
    for (size_t n = 0; n < 3*natom(); n += 3) {
        mass_sum += m_fict_mass[n];
        for (size_t k = 0; k < 3; ++k)
            com[k] += m_fict_mass[n]*m_pos_nmode[n + k];
    }

    for (size_t k = 0; k < 3; ++k)
        com[k] /= mass_sum;

    for (size_t n = 0; n < 3*natom(); n += 3)
        for (size_t k = 0; k < 3; ++k)
            m_pos_nmode[n + k] -= com[k];

    pos_n2c();
}

//----------------------------------------------------------------------------//

double pimd::force(const double* x, double* f)
{
    double Epot = m_potential(nw(), nI(), x, f);

    Epot += m_sphere(nw(), x, f);

    for (size_t n = 0; n < 3*natom(); ++n)
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

    nbead= strtod(argv[1], NULL);

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
	      << "  > nw    = " << sim.nw() << std::endl 
	      << "  > nI    = " << sim.nI() << std::endl 
	      << "  > natom = " << sim.natom() << std::endl ;

    // 2. iterate

    for (size_t n = 0; n < nsteps; ++n) {
        sim.step(dt);
        if (n%nprint == 0) {
            std::cout << n*dt << ' '
                      << sim.invariant() << ' '
                      << sim.Ep() << ' '
                      << sim.Es() << std::endl;
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
        for (size_t w = 0; w < nw(); ++w) {
            print_xyz(ofs, 'O', m_pos_cart + b*ndofs() + 9*w + 0);
            print_xyz(ofs, 'H', m_pos_cart + b*ndofs() + 9*w + 3);
            print_xyz(ofs, 'H', m_pos_cart + b*ndofs() + 9*w + 6);
        }
	print_xyz(ofs, 'I', m_pos_cart + b*ndofs() + 9*nw());
    }
}

////////////////////////////////////////////////////////////////////////////////
