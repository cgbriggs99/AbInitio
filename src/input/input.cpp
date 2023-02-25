#include <cstdio>
#include <stdexcept>
#include "input.hpp"
#include <cctype>
#include "../util/atom.hpp"
#include "../util/molecule.hpp"

using namespace compchem;
using namespace std;

#define BUFF_SIZE 1024

Molecule &compchem::parseXYZ(FILE *fp) {
  /*
   * For this basic implementation, assume that there are no comment
   * lines, and that each line has an element and position.
   */
  char symb[3];
  double x, y, z;
  Molecule *out = new Molecule();

  while(!feof(fp)) {
    char ch = getc(fp);
    if(!isalpha(ch)) {
      continue;
    }
    symb[0] = ch;
    ch = getc(fp);
    symb[1] = (isalpha(ch))? ch: 0;
    symb[2] = 0;
    int got = fscanf(fp, "%lf", &x);
    got += fscanf(fp, "%lf", &y);
    got += fscanf(fp, "%lf", &z);
    if(got != 3 || getZFromSymb(symb) == 0) {
      throw new std::runtime_error("XYZ file not formatted properly.");
    }
    out->addatom(*new Atom(getZFromSymb(symb), x, y, z));
  }
  return *out;
}
  
