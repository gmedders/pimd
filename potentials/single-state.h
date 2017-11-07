#ifndef SINGLE_STATE_H
#define SINGLE_STATE_H

#include <cstdlib>
#include <functional>

#include <armadillo>

#include "explicit-bath.h"
#include "bead-states.h"

namespace pot {

class single_state : public explicit_bath , public bead_states {

public:
    bool init_pot = false;

    double force(size_t, size_t, size_t, const double*, double*);
    virtual double VAA(const double* x, double* f) = 0;
    virtual void set_params(double*) = 0;
    virtual void assert_ndim(int) = 0;

    void set_hopping_params(double*){};

    void print_state_params();

};

} // namespace pot

#endif // SINGLE_STATE_H
