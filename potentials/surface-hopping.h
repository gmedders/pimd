#ifndef SURFACE_HOPPING_H
#define SURFACE_HOPPING_H

#include <cstdlib>

#include <armadillo>

namespace pot {

struct surface_hopping {

    bool init_pot = false;
    bool init_hop = false;
    double Gamma;
    double dt;
    double beta;

    arma::vec fAA;
    arma::vec fBB;

    int active_state;
    int nhops;

    double force(size_t, size_t, const double*, double*);
    virtual double VAA(const double* x, double* f) = 0;
    virtual double VBB(const double* x, double* f) = 0;
    virtual void set_params(double*) = 0;

    void check_allocation(size_t, arma::vec&);

    double fermi_function(const double);
    double hop_probability(const double);
    void hop();

    void set_active_state(const int);
    void set_hopping_params(double*);

};

} // namespace pot

#endif // SURFACE_HOPPING_H
