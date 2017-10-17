#ifndef SURFACE_HOPPING_H
#define SURFACE_HOPPING_H

#include <cstdlib>

#include <armadillo>

#include "explicit-bath.h"

namespace pot {

struct surface_hopping : public explicit_bath {

    bool init_pot = false;
    bool init_hop = false;
    double Gamma;
    double dt;
    double beta;

    arma::ivec state_id;

    int nhops;

    double force(size_t, size_t, size_t, const double*, double*);
    virtual double VAA(const double* x, double* f) = 0;
    virtual double VBB(const double* x, double* f) = 0;
//    virtual double bath_force(const double* x, double* f) = 0;
    virtual void set_params(double*) = 0;

    void check_allocation(size_t, arma::ivec&);

    double fermi_function(const double);
    double hop_probability(const double, int);
    void hop();

    void set_all_bead_states(const int, int);
    void set_individual_bead_states(std::vector<int>&);
    void set_hopping_params(double*);

    double avg_active_state();
    double sum_active_state();

private:
    int m_nbead;

};

} // namespace pot

#endif // SURFACE_HOPPING_H
