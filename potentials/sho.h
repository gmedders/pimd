#ifndef SHO_H
#define SHO_H

#include <cstdlib>

namespace pot {

struct sho {

    double w;
    double m;

    double a;

    double force(const size_t, const size_t,
                 const double* x, double* f);
    double VAA(const double* crd, double* f);

    void set_params(double*);

    bool init = false;

};

} // namespace pot

#endif // SHO_H
