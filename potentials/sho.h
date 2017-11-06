#ifndef SHO_H
#define SHO_H

#include <cstdlib>

#include "single-state.h"

namespace pot {

class sho : public single_state
{

public:
    double VAA(const double* crd, double* f);

    void set_params(double*);

    double get_w(){ return w; };
    double get_m(){ return m; };

    void print_params();

private:
    double w;
    double m;
    double a;
};

} // namespace pot

#endif // SHO_H
