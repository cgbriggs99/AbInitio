#include "../src/input.hpp"
#include <cstdio.h>

using namespace compchem;

int test_input() {
  std::FILE *fp = fopen("water.xyz", "r");

  Molecule *mol = 
