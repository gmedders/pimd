#ifndef NHC_H
#define NHC_H

#include <cstddef>
#include <random>

namespace parts {
namespace nhc {

//
// one degree of freedom assumed (massive thermostating)
//

size_t size(size_t M); // return size (in doubles) of a thermostat

void initialize(size_t M, double *thermo, const double &tau, std::mt19937 &);

// returns velocity scale factor; Ek2 is twice kinetic energy in units of kT
double advance(size_t M, double *thermo, const double &tau, const double &Ek2kT,
               const double &dt);

// in units of kT
double invariant(size_t M, const double *thermo, const double &tau);

} // namespace nhc
} // namespace parts

#endif // NHC_H
