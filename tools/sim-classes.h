#ifndef SIM_CLASSES_H
#define SIM_CLASSES_H

#include "pimd-base.h"
#include "rpmd-base.h"
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
typedef pot::double_well potential_type;
static double omega(2.0e-4); // omega
static double atm_mass(2000); // au
static double bb_x0(21.795);
static double params[] = {omega, atm_mass, bb_x0};
#endif

#if 0
// Anderson-Holstein
typedef pot::ah potential_type;
static double param_w(0.003);
static double atm_mass(2000);
static double param_g(0.02);
static double param_Ed_bar(0.0);
static double params[] = {param_w, atm_mass, param_g, param_Ed_bar};
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
    double avg_cart_pos(void);

    void dump_1D_frame(std::ofstream&);

    potential_type m_potential;

private:
    size_t m_natom;
    size_t m_ndim;
    size_t m_ndofs;
    size_t m_nbead;

    double m_beta;

    double* pos;
};

////////////////////////////////////////////////////////////////////////////////

struct rpmd : public parts::rpmd_base {

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
    double avg_cart_pos(void);

    potential_type m_potential;

private:
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
    double avg_cart_pos(void);

    potential_type m_potential;

private:
    size_t m_natom;
    size_t m_ndim;
    size_t m_ndofs;
};

////////////////////////////////////////////////////////////////////////////////

} // namespace

#endif // SIM_CLASSES_H
