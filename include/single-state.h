#ifndef SINGLE_STATE_H
#define SINGLE_STATE_H

#include <cstdlib>

#include <armadillo>

namespace pot {

struct single_state {

    arma::ivec state_id;

    bool init_pot = false;

    double force(size_t, size_t, size_t, const double*, double*);
    virtual double VAA(const double* x, double* f) = 0;
    virtual void set_params(double*) = 0;

    void set_hopping_params(double*){};

    void set_all_bead_states(const int, int);
    void set_individual_bead_states(std::vector<int>&);

    double avg_active_state(){ return 0; };
    double sum_active_state(){ return 0; };

    void print_state_params();

private:
    int m_nbead;

    void check_allocation(size_t, arma::ivec&);

};

} // namespace pot

#endif // SINGLE_STATE_H
