#ifndef AH_H
#define AH_H

#include <cstdlib>

#include "surface-hopping.h"

namespace pot {

class ah : public surface_hopping
{

public:
    double VAA(const double* x, double* f);
    double VBB(const double* x, double* f);

    void set_params(double*);

    double get_w(){ return w; };
    double get_m(){ return m; };

    void print_params();

private:
    double w;
    double m;
    double g;

    double a;
    double b;
    double Ed;
    double Ed_bar;
};

} // namespace pot

#endif // AH_H
