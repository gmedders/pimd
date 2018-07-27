#ifndef POT_DEFINITION_H
#define POT_DEFINITION_H

#include "ah.h"
#include "anharmonic.h"
#include "double_well.h"
#include "sho.h"

namespace parts {

////////////////////////////////////////////////////////////////////////////////

#if 1
// SHO
typedef pot::sho potential_type;
static double omega(0.001);   // omega
static double atm_mass(2000); // au
static double params[] = {omega, atm_mass};
#endif

#if 0
// DOUBLE WELL
typedef pot::double_well potential_type;
static double omega(0.0009765625); // omega
static double atm_mass(2000); // au
//static double param_g(4.4); // barrier = 3w
static double param_g(3.9); // barrier = 2w
static double dG(-0.003906252);
static double params[] = {omega, atm_mass, param_g, dG};
//
//////static double omega(2.0e-4); // omega
//////static double atm_mass(2000); // au
//////static double param_g(20.6097);
//////static double dG(-0.0038);
//////
//////static double omega(0.001); // omega
//////static double atm_mass(2000); // au
//////static double param_g(3.1);
//////static double dG(-0.004);
//////
//////
//////static double omega(0.006132813); // omega
//////static double atm_mass(2000); // au
//////static double param_g(0.62);
//////static double dG(-0.003906252);
//////
#endif

#if 0
// Anderson-Holstein
typedef pot::ah potential_type;
static double omega(0.003);
static double atm_mass(2000);
static double param_g(0.02);
static double param_Ed_bar(0.1333333);
//static double omega(0.3);
//static double atm_mass(2000);
//static double param_g(0.75);
//static double param_Ed_bar(0.0);
static double params[] = {omega, atm_mass, param_g, param_Ed_bar};
#endif

#if 0
// ANHARMONIC OSCILLATOR
typedef pot::anharmonic potential_type;
static double atm_mass(1); // au
static double param_a(0.0);
static double param_b(0.0);
static double param_c(0.25);
//static double param_a(1.0/2.0);
//static double param_b(1.0/10.0);
//static double param_c(1.0/100.0);
static double params[] = {param_a, param_b, param_c, atm_mass};
#endif

} // namespace parts

#endif // POT_DEFINITION_H
