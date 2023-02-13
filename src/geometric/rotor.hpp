#ifndef __ROTOR_HPP__
#define __ROTOR_HPP__

#include "../util/molecule.hpp"
#include <array>

namespace compchem {

  typedef enum {
    SPHERICAL,
    LINEAR,
    OBLATE,
    PROLATE,
    ASSYMETRIC
  } rotor_t;

  rotor_t comprotor(const Molecule &mol, std::array<double, 3> *out);


}



#endif
