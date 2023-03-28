 
#ifndef __MOLECULE_HPP__
#define __MOLECULE_HPP__

#include "atom.hpp"
#include <vector>

namespace compchem {

  class Molecule {
  private:

    std::vector<Atom *> *atoms;
    double comx, comy, comz;
    int electrons, mult;

    void calc_com();

  public:
    Molecule(const std::vector<Atom *> &atoms, int mult = 1);
    Molecule(const Molecule &copy);
    Molecule();

    virtual ~Molecule();

    Molecule *copy() const;

    int getsize() const;

    const std::vector<Atom *> &getatoms() const;
    const Atom &getatom(int index) const;
    Atom &getatom(int index);

    void addatom(Atom &atom);

    double getcomx() const;
    double getcomy() const;
    double getcomz() const;

    void translate_to_com();

    int getelectrons() const;
    int getmult() const;
  };


  double nuclear_repulsion(const Molecule &mol);
}

#endif
