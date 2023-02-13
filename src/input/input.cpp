#include <cstdio>
#include "input.hpp"
#include "../util/atom.hpp"
#include "../util/molecule.hpp"

using namespace compchem;

#define BUFF_SIZE 1024

Molecule &compchem::parseXYZ(std::FILE *fp) {
  /*
   * For this basic implementation, assume that there are no comment
   * lines, and that each line has an element and position.
   */
  char symb[3];
  double x, y, z;
  Molecule *out = new Molecule();

  while(!std::feof(fp)) {
    fscanf(fp, "%s %d %d %d\n", symb, &x, &y, &z);

    out->addatom(*new Atom(getZFromSymb(symb), x, y, z));
  }
  return *out;
}
  
