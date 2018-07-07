#ifndef GCMC_H
#define GCMC_H

#include <memory>

#include "sim-classes.h"

//
// units are au
//

namespace parts {

////////////////////////////////////////////////////////////////////////////////

class gcmc {

public:
  gcmc() : m_md_ensemble(new rpmd){};

  void set_chemical_potential(double);
  // void set_md_ensemble(std::string&);
  bool calc_insertion_probability();

  void set_up(const size_t, const size_t, const size_t, const double,
              const double, double *, double *);

  void step(double, double);
  void dump(std::ostream &);

  // typedef ensemble_type m_md_ensemble;
  std::unique_ptr<rpmd> m_md_ensemble;

private:
  double m_chemical_potential = 0;
};

////////////////////////////////////////////////////////////////////////////////

} // namespace parts

#endif // GCMC_H
