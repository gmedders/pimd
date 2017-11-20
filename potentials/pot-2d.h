#ifndef POT_2D_H
#define POT_2D_H

#include <cstdlib>

#include "single-state.h"

namespace pot {

class pot_2d : public single_state
{

public:
    //pot_2d() {std::cerr << "pot_2d::pot_2d()" << std::endl; };
    //~pot_2d() {std::cerr << "pot_2d::~pot_2d()" << std::endl; };
    double VAA(const double* crd, double* f);
    double VBB(const double* crd, double* f);

    void assert_ndim(int);

    void set_params(double*);

    double get_w(){ return w; };
    double get_m(){ return mx; };

    void print_params();

private:
    double A0, B1;
    double x0, x1;
    double C0, C1;
    double w, mx, z0;
};

} // namespace pot

#endif // POT_2D_H
