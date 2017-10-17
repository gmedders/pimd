#ifndef EXPLICIT_BATH_H
#define EXPLICIT_BATH_H

#include <cstdlib>

#include <armadillo>

namespace pot {

struct explicit_bath {

    explicit_bath();
    ~explicit_bat();

    void init(size_t ndim, size_t nBathModes,
              double gamma, double cutoff);
    double bath_force(const double*, double*);

private:
    int m_nBathModes;

    double* c_iMode;
    double* omega_iMode; 

};

} // namespace pot

#endif // EXPLICIT_BATH_H
