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

//
// units are au
//

namespace parts {

////////////////////////////////////////////////////////////////////////////////

#if 0
// SHO
typedef pot::sho potential_type;
static double omega(0.2); // omega
static double atm_mass(2000); // au
static double params[] = {omega, atm_mass};
#endif

#if 1
// DOUBLE WELL
//typedef pot::double_well potential_type;
//static double omega(0.0009765625); // omega
//static double atm_mass(2000); // au
//static double param_g(3.9);
//static double dG(-0.003906252);
//static double params[] = {omega, atm_mass, param_g, dG};
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

#if 1
// Anderson-Holstein
typedef pot::ah potential_type;
//static double omega(0.003);
//static double atm_mass(2000);
//static double param_g(0.02);
//static double param_Ed_bar(0.0);
static double omega(0.3);
static double atm_mass(2000);
static double param_g(0.75);
static double param_Ed_bar(0.0);
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
static double params[] = {param_a, param_b, param_c};
#endif

////////////////////////////////////////////////////////////////////////////////

struct pimd : public parts::pimd_base {

    ~pimd();

    void set_up(const size_t, const size_t, const size_t, const double);
    double force(const size_t, const size_t, const double*, double*);

    inline double Espring() const { return m_Espring; }
    inline double Ep() const { return m_Epot_sum; }
    inline double Ek() const { return m_Ekin_fict; }
    inline double temp_kT() const { return m_temp_kT; }
    double avg_cart_pos() { calc_pos_stats(); return m_avg_cart_pos; };
    double L1_cart_pos() const { return m_L1_cart_pos; };
    double L2_cart_pos() const { return m_L2_cart_pos; };
    double Linf_cart_pos() const { return m_Linf_cart_pos; };
    void calc_pos_stats(void);

    void dump_1D_frame(std::ofstream&);

    potential_type m_potential;

private:
    size_t m_natom;
    size_t m_ndim;
    size_t m_ndofs;
    size_t m_nbead;

    double m_beta;

    double m_avg_cart_pos;
    double m_L1_cart_pos;
    double m_L2_cart_pos;
    double m_Linf_cart_pos;

    double* pos;
};

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
    double force(const size_t, const size_t, const double*, double*);

    inline double Espring() const { return m_Espring; }
    inline double Ep() const { return m_Epot_sum; }
    inline double Ek() const { return m_Ekin; }
    inline double temp_kT() const { return m_temp_kT; }
    double avg_cart_pos() { calc_pos_stats(); return m_avg_cart_pos; };
    //double avg_cart_pos() const { return m_avg_cart_pos; };
    double L1_cart_pos() const { return m_L1_cart_pos; };
    double L2_cart_pos() const { return m_L2_cart_pos; };
    double Linf_cart_pos() const { return m_Linf_cart_pos; };
    void calc_pos_stats(void);

    void dump_1D_frame(std::ofstream&);

    potential_type m_potential;

private:
    //double gamma = 5.0*2.0*omega;
    double gamma = 2.0*omega;

    double m_avg_cart_pos;
    double m_L1_cart_pos;
    double m_L2_cart_pos;
    double m_Linf_cart_pos;

    size_t m_natom;
    size_t m_ndim;
    size_t m_ndofs;
    size_t m_nbead;
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
    double force(const size_t, const size_t, const double*, double*);

    inline double Espring() const { return 0.0; }
    inline double Ep() const { return m_Epot; }
    inline double Ek() const { return m_Ekin; }
    inline double temp_kT() const { return m_temp_kT; }
    double avg_cart_pos() { calc_pos_stats(); return m_avg_cart_pos; };
    //double avg_cart_pos() const { return m_avg_cart_pos; };
    double L1_cart_pos() const { return m_L1_cart_pos; };
    double L2_cart_pos() const { return m_L2_cart_pos; };
    double Linf_cart_pos() const { return m_Linf_cart_pos; };
    void calc_pos_stats(void);

    potential_type m_potential;

private:
    double m_avg_cart_pos;
    double m_L1_cart_pos;
    double m_L2_cart_pos;
    double m_Linf_cart_pos;

    size_t m_natom;
    size_t m_ndim;
    size_t m_ndofs;
};

////////////////////////////////////////////////////////////////////////////////

} // namespace

#endif // SIM_CLASSES_H
