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

    int active_state;
    void set_active_state(const int);
    void set_hopping_params(double*);

};

} // namespace pot

#endif // SHO_H
