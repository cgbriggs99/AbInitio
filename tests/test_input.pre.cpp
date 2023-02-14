#include "${CMAKE_SOURCE_DIR}/src/input/input.hpp"
#include "${CMAKE_SOURCE_DIR}/src/util/molecule.hpp"
#include <cstdio>
#include <errno.h>

using namespace compchem;

int test_input() {
  errno = 0;
  std::FILE *fp = std::fopen("${CMAKE_SOURCE_DIR}/tests/water.xyz", "r");
  if(fp == nullptr || errno != 0) {
    std::perror("Could not open ${CMAKE_SOURCE_DIR}/tests/water.xyz");
    return 1;
  }

  Molecule *mol = &compchem::parseXYZ(fp);
  
  std::fclose(fp);

  int fails = 0;
  if(mol->getsize() != 3) {
    fails++;
  }
  delete mol;
  return fails;
}

int main(void) {
  int fails = 0;
  fails += test_input();

  return fails;
}
