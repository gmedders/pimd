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

};

} // namespace pot

#endif // ANHARMONIC_H
