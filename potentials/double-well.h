#ifndef DOUBLE_WELL_H
#define DOUBLE_WELL_H

#include <cstdlib>

#include <armadillo>

#include "mt19937.h"

namespace pot {

struct double_well {

    bool init_pot = false;
    double w;
    double m;
    double g;

    bool init_hop = false;
    double Gamma;
    double dt;
    double beta;

    double a;

    arma::vec fAA;
    arma::vec fBB;

    int active_state;
    int nhops;

    double force(size_t, size_t, const double*, double*);
    double VAA(const double* x, double* f);
    double VBB(const double* x, double* f);

    void check_allocation(size_t, arma::vec&);

    double fermi_function(const double);
    double hop_probability(const double);
    void hop();

    void set_params(double*);
    void set_hopping_params(double*);

    parts::mt19937 prng;

};

} // namespace pot

#endif // DOUBLE_WELL_H
