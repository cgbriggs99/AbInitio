#ifndef __PSI4_BS_READER_HPP__
#define __PSI4_BS_READER_HPP__

#include "../basis_sets/basis_set.hpp"
#include <vector>
#include <cstdio>

namespace compchem {

  std::vector<GaussianOrbital> *readPsi4file(std::FILE *fp, int Z, int charge = 0);

}


#endif
