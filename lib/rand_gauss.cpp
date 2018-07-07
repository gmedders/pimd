#include <math.h>
#include <stdlib.h>

namespace parts {

double randn(double mu, double sigma) {
  double U1, U2, W, mult;
  static double X1, X2;
  static bool first_call = true;

  if (!first_call) {
    first_call = true;
    return (mu + sigma * (double)X2);
  }

  do {
    U1 = -1 + ((double)rand() / RAND_MAX) * 2;
    U2 = -1 + ((double)rand() / RAND_MAX) * 2;
    W = pow(U1, 2) + pow(U2, 2);
  } while (W >= 1 || W == 0);

  mult = sqrt((-2 * log(W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;

  first_call = false;

  return (mu + sigma * (double)X1);
}

} // namespace parts
