
#ifndef __ATOM_HPP__
#define __ATOM_HPP__

#include "../basis_sets/basis_set.hpp"
#include "copyable.hpp"
#include <string>
#include <vector>

namespace compchem {

  class Atom : Copyable {
  private :
    int Z;
    int charge;
    double x, y, z, mass;
    int norbitals;
    std::vector<BasisOrbital *> orbitals;
  public :
    Atom();
    Atom(int Z, double x, double y, double z);
    Atom(int Z, double mass, double x, double y, double z);
    Atom(int Z, double x, double y, double z, int norbs, const BasisOrbital *passorbs);
    Atom(int Z, double mass, double x, double y, double z, int norbs, const BasisOrbital *passorbs);
    Atom(int Z, int charge, double x, double y, double z);
    Atom(int Z, int charge, double mass, double x, double y, double z);
    Atom(int Z, int charge, double x, double y, double z, int norbs, const BasisOrbital *passorbs);
    Atom(int Z, int charge, double mass, double x, double y, double z, int norbs, const BasisOrbital *passorbs);

    // Copy constructor.
    Atom(const Atom &copy);

    virtual ~Atom();

    double getx() const;
    double gety() const;
    double getz() const;

    double getmass() const;
    int getZ() const;
    int getcharge() const;
    int getnorbitals() const;
    const std::vector<BasisOrbital *> &getorbitals() const;
    const BasisOrbital &getorbital(int index) const;

    Atom &copy() override;
  };

  int getZFromSymb(const std::string &symb);
  double getAbundantMass(int Z);
}


#endif
