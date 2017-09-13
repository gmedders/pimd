#ifndef SIM_CLASSES_H
#define SIM_CLASSES_H

#include "pimd-base.h"
#include "rpmd-base.h"

#include "sho.h"
#include "anharmonic.h"

//
// units are au
//

namespace parts {

////////////////////////////////////////////////////////////////////////////////

//typedef pot::sho potential_type;
//static double atm_mass(2000); // au
//static double tomega(0.2); // omega
//static double params[] = {omega, atm_mass};

typedef pot::anharmonic potential_type;
static double atm_mass(1); // au
static double param_a(0.0);
static double param_b(0.0);
static double param_c(0.25);
//static double param_a(1.0/2.0);
//static double param_b(1.0/10.0);
//static double param_c(1.0/100.0);
static double params[] = {param_a, param_b, param_c};

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

private:
    size_t m_natom;
    size_t m_ndim;
    size_t m_ndofs;
    size_t m_nbead;

    double m_beta;

    double* pos;
    potential_type m_potential;
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

private:
    size_t m_natom;
    size_t m_ndim;
    size_t m_ndofs;
    size_t m_nbead;
    potential_type m_potential;
};

////////////////////////////////////////////////////////////////////////////////

} // namespace

#endif // SIM_CLASSES_H
