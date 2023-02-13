#ifndef __COMPCHEM_INPUT_HPP__
#define __COMPCHEM_INPUT_HPP__

#include "../util/molecule.hpp"
#include "../util/atom.hpp"
#include <cstdio>
#include <vector>

namespace compchem {

  Molecule &parseXYZ(std::FILE *fp);

}


#endif
