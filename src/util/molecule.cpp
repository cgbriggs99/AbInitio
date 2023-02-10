 
#include "molecule.hpp"
#include <vector>

using namespace compchem;

Molecule::Molecule(const std::vector<Atom *> &atoms) {
  this->atoms = std::vector<Atom *>(atoms.size());

  for(Atom *atom : atoms) {
    this->atoms.push_back(atom->copy());
  }
}

Molecule::Molecule(const Molecule &copy) {
  this->atoms = std::vector<Atom *>(copy.atoms.size());

  for(Atom *atom : copy.atoms) {
    this->atoms.push_back(atom->copy());
  }
}

Molecule::Molecule() {
  this->atoms = std::vector<Atom *>(0);
}

Molecule::~Molecule() {
  for(Atom *atom : this->atoms) {
    delete atom;
  }
}

Molecule &Molecule::copy() {
  return *new Molecule(*this);
}

int Molecule::getsize() {
  return this->atoms.size();
}

const std::vector<Atom *> &Molecule::getatoms() const {
  return this->atoms;
}

const Atom &getatom(int index) const {
  return *this->atoms.at(index);
}
