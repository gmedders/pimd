#ifndef SIM_CLASSES_H
#define SIM_CLASSES_H

#include "pimd-base.h"

#include "rpmd-base.h"
#include "rpmd-nhc.h"
#include "rpmd-pile.h"

#include "vv-base.h"

#include "sho.h"
#include "anharmonic.h"
#include "double-well.h"
#include "ah.h"
#include "pot-2d.h"

//
// units are au
//

namespace parts {

////////////////////////////////////////////////////////////////////////////////

#if 0
// SHO
typedef pot::sho potential_type;
static double omega(0.001); // omega
static double atm_mass(2000); // au
static double params[] = {omega, atm_mass};
#endif

#if 1
// pot_2d
typedef pot::pot_2d potential_type;
static double A0 = 0.011;
static double B1 = 0.2;
static double x0 = 0.0;
static double x1 = 0.18899;
static double C0 = 0.64;
static double C1 = 0.67;
static double w  = 8.0E-3;
static double mx = 14000;
static double z0 = -3.5;

static double atm_mass = 55000;

static double params[] = {
A0,
B1,
x0,
x1,
C0,
C1,
w ,
mx,
z0};
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

////////////////////////////////////////////////////////////////////////////////

//struct rpmd : public parts::rpmd_base {
struct rpmd : public parts::rpmd_pile {
//struct rpmd : public parts::rpmd_nhc {

    void set_up_new_init_cond(const size_t, const size_t, const size_t,
                              const double, const double);
    void set_up_new_init_cond(const size_t, const size_t, const size_t,
                              const double, const double, double*);
    void set_up(const size_t, const size_t, const size_t,
                const double, const double,
                double*, double*);
    double force(size_t, size_t, size_t, const double*, double*);

    inline double Espring() const { return m_Espring; }
    inline double Ep() const { return m_Epot_sum; }
    inline double Ek() const { return m_Ekin; }
    inline double temp_kT() const { return m_temp_kT; }
    inline double temp_kT_centroid() const { return m_temp_kT_centroid; }
    inline double temp_kT_higherNM() const { return m_temp_kT_higherNM; }
    double avg_cart_pos() { calc_pos_stats(); return m_avg_cart_pos; };
    //double avg_cart_pos() const { return m_avg_cart_pos; };
    double L1_cart_pos() const { return m_L1_cart_pos; };
    double L2_cart_pos() const { return m_L2_cart_pos; };
    double Linf_cart_pos() const { return m_Linf_cart_pos; };
    void calc_pos_stats(void);

    void set_gammaTh(const double&, double);

    void dump_1D_frame(std::ofstream&);
    void print_params();

    potential_type m_potential;

    double m_gamma = 0.0;

private:
    double m_avg_cart_pos;
    double m_L1_cart_pos;
    double m_L2_cart_pos;
    double m_Linf_cart_pos;
};

////////////////////////////////////////////////////////////////////////////////

struct vv : public parts::vv_base {

    void set_up_new_init_cond(const size_t, const size_t, const size_t,
                              const double, const double);
    void set_up_new_init_cond(const size_t, const size_t, const size_t,
                              const double, const double, double*);
    void set_up(const size_t, const size_t, const size_t,
                const double, const double,
                double*, double*);
    double force(size_t, size_t, size_t, const double*, double*);

    inline double Espring() const { return 0.0; }
    inline double Ep() const { return m_Epot; }
    inline double Ek() const { return m_Ekin; }
    inline double temp_kT() const { return m_temp_kT; }

    double avg_cart_pos() { calc_pos_stats(); return m_avg_cart_pos; };
    double L1_cart_pos() const { return 0.0; };
    double L2_cart_pos() const { return 0.0; };
    double Linf_cart_pos() const { return 0.0; };
    void calc_pos_stats(void);

    void print_params();

    potential_type m_potential;

private:
    double m_avg_cart_pos;
};

////////////////////////////////////////////////////////////////////////////////

} // namespace

#endif // SIM_CLASSES_H
