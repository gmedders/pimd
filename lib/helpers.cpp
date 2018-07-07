#include "helpers.h"

#include <iostream>
#include <sstream>

namespace parts {

//----------------------------------------------------------------------------//

int parse_to_int(char *argv) {
  int val;
  std::istringstream iss(argv);
  iss >> val;
  if (!iss || !iss.eof()) {
    std::cerr << "could not convert '" << argv << "' to int" << std::endl;
    return EXIT_FAILURE;
  }
  return val;
}

//----------------------------------------------------------------------------//

double parse_to_double(char *argv) {
  double val;
  std::istringstream iss(argv);
  iss >> val;
  if (!iss || !iss.eof()) {
    std::cerr << "could not convert '" << argv << "' to double" << std::endl;
    return EXIT_FAILURE;
  }
  return val;
}

//----------------------------------------------------------------------------//

} // namespace parts
