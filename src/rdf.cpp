#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <cassert>
#include <cmath>
#include <cstdlib>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "sim-classes.h"

////////////////////////////////////////////////////////////////////////////////

namespace {

const double r_min = -50;
const double r_max = 50;
const double dr = 0.01;

void check_parsing(std::istringstream &iss, size_t lineno) {
  if (iss.fail()) {
    std::ostringstream oss;
    oss << "failed to parse line " << lineno
        << " of the input file:" << std::endl
        << iss.str() << std::endl;
    throw std::runtime_error(oss.str());
  }
}

} // namespace

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
  // 1. load the coordinates

  std::cout.setf(std::ios_base::showpoint);
  std::cout.precision(9);

  if (argc != 2) {
    std::cerr << "usage: rdf input_file" << std::endl;
    return EXIT_FAILURE;
  }

  size_t ndim = 1;
  size_t natom = 1;

  // 2. iterate
  std::string filename(argv[1]);
  std::ifstream ifs(filename.c_str());

  size_t lineno(0);

  // Initialize arrays
  double nsamples(0);
  std::vector<double> rdf;
  std::vector<double> pos;

  int nbin = int((r_max - r_min) / dr);

  for (size_t i = 0; i < nbin; ++i) {
    const double r = dr * (i + 0.5) + r_min;

    rdf.push_back(0.0);
    pos.push_back(r);
  }

  while (!ifs.eof()) {

    // Read this frame of the input file
    // Fist line: NBead NDof Beta
    std::string line;
    std::getline(ifs, line);
    ++lineno;

    if (ifs.eof())
      break;

    size_t nbead;
    size_t ndof;
    double beta;
    size_t init_active_state;
    std::istringstream iss(line);
    iss >> nbead >> ndof >> beta >> init_active_state;
    check_parsing(iss, lineno);

    assert(nbead > 0);
    assert(ndof == ndim * natom);
    assert(beta > 0);

    for (size_t n = 0; n < nbead; ++n) {
      std::string line;
      std::getline(ifs, line);
      ++lineno;
      std::istringstream iss(line);

      for (size_t i = 0; i < ndof; ++i) {
        double q;
        double v;

        iss >> q >> v;
        check_parsing(iss, lineno);

        int binid = int((q - r_min) / dr);
        //       std::cout << q << ' ' << binid << ' ' << rdf.size() <<
        //       std::endl;

        assert(binid < rdf.size());
        ++rdf[binid];
        ++nsamples;
      }
    }
  }

  // Finally, print the results
  std::cout << std::scientific;
  std::cout.precision(10);
  for (size_t i = 0; i < pos.size(); ++i) {
    std::cout << std::setw(20) << pos[i] << std::setw(20) << rdf[i] / nsamples
              << std::endl;
  }

  ifs.close();

  return EXIT_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
