#ifndef BEAD_STATES_H
#define BEAD_STATES_H

#include <cstdlib>

#include <armadillo>

namespace pot {

class bead_states {

public:
  // bead_states() {std::cerr << "bead_states::bead_states()" << std::endl; };
  //~bead_states() {std::cerr << "bead_states::~bead_states()" << std::endl; };
  arma::ivec state_id;

  void set_all_bead_states(const int, int);
  void set_individual_bead_states(std::vector<int> &);

  double avg_active_state();
  double sum_active_state();

  void check_allocation(size_t, arma::ivec &);

private:
  int m_nbead;
};

} // namespace pot

#endif // BEAD_STATES_H
