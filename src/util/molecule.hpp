 
#ifndef __MOLECULE_HPP__
#define __MOLECULE_HPP__

#include "copyable.hpp"
#include "atom.hpp"
#include <vector>

namespace compchem {

  class Molecule : Copyable {
  private:

    std::vector<Atom *> atoms;

  public:
    Molecule(const std::vector<Atom *> &atoms);
    Molecule(const Molecule &copy);
    Molecule();

    virtual ~Molecule();

    Molecule &copy() override;

    int getsize() const;

    const std::vector<Atom *> &getatoms() const;
    const Atom &getatom(int index) const;
  };
    
}


#endif
