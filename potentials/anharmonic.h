#ifndef ANHARMONIC_H
#define ANHARMONIC_H

#include <cstdlib>

namespace pot {

struct anharmonic {

    double a;
    double b;
    double c;

    double force(const size_t, const size_t,
                 const double* x, double* f);
    double VAA(const double* crd, double* f);

    void set_params(double*);

    bool init = false;

    int active_state;
    void set_all_bead_states(const int);
    void set_hopping_params(double*);

};

} // namespace pot

#endif // ANHARMONIC_H
