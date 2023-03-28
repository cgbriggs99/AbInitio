#ifndef WAVEFUNCTION_HPP
#define WAVEFUNCTION_HPP

#include "../util/molecule.hpp"

namespace compchem {

class Wavefunction {
protected:
  int electrons, multiplicity;
public:
  Wavefunction(int electrons) : electrons(electrons),
				multiplicity(1)
  {;}
  Wavefunction(int electrons, int multiplicity) : electrons(electrons),
						  multiplicity(multiplicity)
  {;}

  virtual ~Wavefunction() = default;
  
  virtual const double *getoverlap(int *dim = nullptr) const = 0;
  virtual const double *getkinetic(int *dim = nullptr) const = 0;
  virtual const double *getpotential(int *dim = nullptr) const = 0;

  virtual const double *getfock(int *dim) const = 0;

  virtual int getnorbs() const = 0;

  virtual int getelectrons() const;
  virtual int getmultiplicity() const;

};

}



#endif
