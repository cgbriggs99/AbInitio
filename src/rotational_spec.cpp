#include "${CMAKE_SOURCE_DIR}/src/scf/rhf.hpp"
#include "${CMAKE_SOURCE_DIR}/src/input/psi4_bs_reader.hpp"
#include "${CMAKE_SOURCE_DIR}/src/basis_sets/sto-ng/stong.hpp"
#include "${CMAKE_SOURCE_DIR}/src/basis_sets/basis_set.hpp"
#include "${CMAKE_SOURCE_DIR}/tests/test.h"
#include "${CMAKE_SOURCE_DIR}/src/integrals/integrals.hpp"
#include "${CMAKE_SOURCE_DIR}/src/input/input.hpp"
#include "${CMAKE_SOURCE_DIR}/src/opts/default_options.hpp"
#include <cstdio>
#include <string>
#include <errno.h>

using namespace compchem;
using namespace std;

void print_peaks(

int main(int argc, char **argv) {
  double temp = 298.15;

  FILE *geo = nullptr, *basis = nullptr;;

  for(int i = 0; i < argc; i++) {
    if(strncmp(argv[i], "-h", 3)) {
      printf("-h: Help\n"
	     "-t: Temperature in Kelvin.\n"
	     "-geo: File with atomic coordinates.\n"
	     "-bas: File for the basis set.\n");
      return;
    } else if(strncmp(argv[i], "-t", 3)) {
      if(i == argc - 1) {
	perror("Missing temperature.");
	exit(1);
      }
      sscanf(argv[i + 1], "%lf", &temp);
      i++;
    } else if(strncmp(argv[i], "-geo", 5)) {
      if(i == argc - 1) {
	perror("Missing geometry file.");
	exit(1);
      }

      geo = fopen(argv[i + 1], "r");
      i++;
    } else if(strncmp(argv[i], "-bas", 5)) {
      if(i == argc - 1) {
	perror("Missing basis file.");
	exit(1);
      }

      basis = fopen(argv[i + 1], "r");
      i++;
    } else {
      perror("Unrecognized argument.");
      exit(1);
    }

    if(geo == nullptr) {
      geo = fopen("${CMAKE_SOURCE_DIR}/tests/water2.xyz", "r");
    }
    if(basis == nullptr) {
      basis = fopen("${CMAKE_SOURCE_DIR}/tests/sto-3g.gbs", "r");
    }

    
