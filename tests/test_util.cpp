
#include "../src/util/atom.hpp"
#include "../src/util/molecule.hpp"
#include <vector>

using namespace compchem;

int test_atom() {
  int fails = 0;
  if(getZFromSymb("Au") != 79) {
    fails++;
  }
  if(getZFromSymb("D") != 1) {
    fails++;
  }
  if(getZFromSymb("T") != 1) {
    fails++;
  }
  if(getZFromSymb("A") != 0) {
    fails++;
  }
  return fails;
}

int test_molecule() {
  Atom *H1 = new Atom(1, 0, 1.42396298, 1.11787308),
    *H2 = new Atom(1, 0, -1.42396298, 1.11787308),
    *O = new Atom(8, 0, 0, 0);

  int fails = 0;
  if(H1->getmass() <= 1) {
    fails++;
  }
  if(H2->getmass() <= 1) {
    fails++;
  }
  if(O->getmass() <= 15) {
    fails++;
  }
  std::vector<Atom *> *atoms = new std::vector<Atom *>({H1, H2, O});
  Molecule *mol = new Molecule(*atoms);

  if(mol->getcomx() != 0) {
    fails++;
  }
  if(mol->getcomy() != 0) {
    fails++;
  }
  if(mol->getcomz() <= 0) {
    fails++;
  }
  mol->translate_to_com();
  if(mol->getcomx() != 0) {
    fails++;
  }
  if(mol->getcomy() != 0) {
    fails++;
  }
  if(mol->getcomz() != 0) {
    fails++;
  }
  delete mol;
  delete atoms;
  delete H1;
  delete H2;
  delete O;
  return fails;
}


int main(void) {
  int fails = 0;

  fails += test_atom();
  fails += test_molecule();


  return fails;

}
