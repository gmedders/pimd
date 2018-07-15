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

#include <armadillo>

#include "helpers.h"
#include "sim_classes.h"

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

  if (argc < 3) {
    std::cerr << "usage: rdf Nframes Nrdfs traj.dat traj.dat ..." << std::endl;
    return EXIT_FAILURE;
  }

  int Nframes = parts::parse_to_int(argv[1]);
  int Nrdfs = parts::parse_to_int(argv[2]);
  assert(Nframes > Nrdfs);

  int NframesPerRDF = int(Nframes / Nrdfs);

  size_t ndim = 1;
  size_t natom = 1;

  // Initialize arrays
  int nbin = int((r_max - r_min) / dr);

  arma::vec pos(nbin, arma::fill::zeros);
  arma::vec nsamples(Nrdfs, arma::fill::zeros);
  arma::mat rdfs(Nrdfs, nbin, arma::fill::zeros);

  for (size_t i = 0; i < nbin; ++i) {
    const double r = dr * (i + 0.5) + r_min;

    pos(i) = r;
  }

  // Now read each file
  for (int ifile = 3; ifile < argc; ++ifile) {
    std::string filename(argv[ifile]);
    std::ifstream ifs(filename.c_str());

    size_t lineno(0);
    size_t iframe(0);
    size_t which_rdf(0);

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
      std::istringstream iss(line);
      iss >> nbead >> ndof >> beta;
      check_parsing(iss, lineno);

      assert(nbead > 0);
      assert(ndof == ndim * natom);
      assert(beta > 0);

      for (size_t n = 0; n < nbead; ++n) {
        std::string line;
        std::getline(ifs, line);
        ++lineno;
        std::istringstream iss(line);

        int this_state;
        iss >> this_state;

        for (size_t i = 0; i < ndof; ++i) {
          double q;
          double v;

          iss >> q >> v;
          check_parsing(iss, lineno);

          ////if(q > r_min && q < r_max){

          int binid = int((q - r_min) / dr);

          rdfs(which_rdf, binid) += 1.0;
          nsamples(which_rdf) += 1.0;
          //}
        }
      }
      ++iframe;
      // if(iframe > 1)
      //    break;
      which_rdf = (int)std::floor(iframe / NframesPerRDF);
    }
    ifs.close();
  }

  // Finally, print the results
  std::cout << std::scientific;
  std::cout.precision(10);
  for (size_t i = 0; i < nbin; ++i) {
    std::cout << std::setw(20) << pos(i);
    for (size_t irdf = 0; irdf < Nrdfs; ++irdf) {
      std::cout << std::setw(20) << rdfs(irdf, i) / nsamples(irdf);
    }
    std::cout << std::endl;
  }

  return EXIT_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
