
#include "../src/util/atom.hpp"
#include "../src/util/molecule.hpp"
#include <vector>

int test_atom() {
  Atom &atom1 = *new Atom();
  Atom &atom2 = *new Atom(1, 0, 1, 2);
  Atom &atom3 = *new Atom(2, 4.0, 0, 1, 2);
  Atom &atom4 = *new Atom(3, 0, 1, 2, 

}



int main(void) {
  int fails = 0;



  return fails;

}
