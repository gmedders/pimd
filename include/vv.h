#ifndef VV_H
#define VV_H

#include <armadillo>

#include "classical.h"

namespace parts {

struct vv : public classical {
  void step(const double &);
  double invariant() const;
};

} // namespace parts

#endif // VV_H
