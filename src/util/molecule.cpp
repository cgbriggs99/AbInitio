 
#include "molecule.hpp"
#include <vector>
#include <iostream>

using namespace compchem;

Molecule::Molecule(const std::vector<Atom *> &atoms) {
  this->atoms = new std::vector<Atom *>(0);

  for(Atom *atom : atoms) {
    Atom *copy = new Atom(*atom);
    this->atoms->push_back(copy);
  }
  this->comx = 0;
  this->comy = 0;
  this->comz = 0;
  this->calc_com();
}

Molecule::Molecule(const Molecule &copy) {
  this->atoms = new std::vector<Atom *>(0);

  for(Atom *atom : *(copy.atoms)) {
    this->atoms->push_back(new Atom(*atom));
  }
  this->comx = copy.comx;
  this->comy = copy.comy;
  this->comz = copy.comz;
}

Molecule::Molecule() {
  this->atoms = new std::vector<Atom *>(0);
  this->comx = 0;
  this->comy = 0;
  this->comz = 0;
}

Molecule::~Molecule() {
  for(Atom *atom : *this->atoms) {
    delete atom;
  }
  delete this->atoms;
}

Molecule *Molecule::copy() const {
  return new Molecule(*this);
}

int Molecule::getsize() const {
  return this->atoms->size();
}

const std::vector<Atom *> &Molecule::getatoms() const {
  return *(this->atoms);
}

const Atom &Molecule::getatom(int index) const {
  return *(this->atoms->at(index));
}

Atom &Molecule::getatom(int index) {
  return *(this->atoms->at(index));
}

void Molecule::addatom(Atom &atom) {
  this->atoms->push_back(&atom);
  this->calc_com();
}

void Molecule::calc_com() {
  double norm = 0, x = 0, y = 0, z = 0;
  if(this->getsize() == 0) {
    this->comx = 0;
    this->comy = 0;
    this->comz = 0;
    return;
  }
  for(Atom *atom : *(this->atoms)) {
    norm += atom->getmass();
    x += atom->getmass() * atom->getx();
    y += atom->getmass() * atom->gety();
    z += atom->getmass() * atom->getz();
  }
  this->comx = x / norm;
  this->comy = y / norm;
  this->comz = z / norm;
}

void Molecule::translate_to_com() {
  for(int i = 0; i < this->getsize(); i++) {
    Atom &atom = this->getatom(i);
    atom.setx(atom.getx() - this->comx);
    atom.sety(atom.gety() - this->comy);
    atom.setz(atom.getz() - this->comz);
  }
  this->comx = 0;
  this->comy = 0;
  this->comz = 0;
}

double Molecule::getcomx() const {
  return this->comx;
}

double Molecule::getcomy() const {
  return this->comy;
}

double Molecule::getcomz() const {
  return this->comz;
}
