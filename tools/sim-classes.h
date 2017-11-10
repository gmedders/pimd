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

#include "pot-definition.h"

//
// units are au
//

namespace parts {

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
