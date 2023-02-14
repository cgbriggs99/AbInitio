#include "${CMAKE_SOURCE_DIR}/src/input/input.hpp"
#include "${CMAKE_SOURCE_DIR}/src/util/molecule.hpp"
#include "${CMAKE_SOURCE_DIR}/src/geometric/rotor.hpp"
#include <cstdio>
#include <array>

using namespace compchem;
using namespace std;

int test_rotor() {
  FILE *fp = fopen("${CMAKE_SOURCE_DIR}/tests/water.xyz", "r");
  if(fp == nullptr) {
    perror("Could not open ${CMAKE_SOURCE_DIR}/tests/water.xyz");
    return 1;
  }

  Molecule *mol = &compchem::parseXYZ(fp);
  array<double, 3> eigs;
  

  fclose(fp);

  int fails = 0;
  if(comprotor(*mol, &eigs) != ASSYMETRIC) {
    fails++;
  }

  for(int i = 0; i < 3; i++) {
    if(eigs[i] <= 0) {
      fails++;
    }
  }
  delete mol;
  return fails;
}

int main(void) {
  int fails = 0;

  fails += test_rotor();

  return fails;
}
