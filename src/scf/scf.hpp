#ifndef SCF_HPP
#define SCF_HPP

#include "../opts/options.hpp"
#include "../wavefunction/wavefunction.hpp"
#include "../util/tei_array.hpp"
#include "../energy/energy.hpp"
#include "../util/molecule.hpp"

namespace compchem {

class SCFWfn : public Wavefunction {
protected:
  int dim;
  const double *S, *T, *V;
  const TEIArray *tei;

  double *Ca, *Cb, *Da, *Db, *Fa, *Fb;
  double energy;
  double *es;
  int orbs;
public:
  
  SCFWfn(int electrons, const double *S, const double *T,
	 const double *V,
	 const TEIArray *tei, int orbs, double *Ca = nullptr,
	 double *Cb = nullptr,
	 double *Da = nullptr,
	 double *Db = nullptr,
	 double *Fa = nullptr,
	 double *Fb = nullptr,
	 double *es = nullptr,
	 int multiplicity = 1) : S(S), T(T), V(V),
				 tei(tei),
				 Ca(Ca),
				 Cb(Cb),
				 Da(Da),
				 Db(Db),
				 Fa(Fa),
				 Fb(Fb),
				 energy(0),
				 es(es),
				 Wavefunction(electrons, multiplicity),
				 orbs(orbs) {;}
  virtual ~SCFWfn();

  virtual const double *getoverlap(int *dim = nullptr) const override;
  virtual const double *getkinetic(int *dim = nullptr) const override;
  virtual const double *getpotential(int *dim = nullptr) const override;

  virtual const TEIArray *gettei() const;
  virtual const double *getcoefa(int *dim = nullptr) const;
  virtual const double *getcoefb(int *dim = nullptr) const;
  virtual const double *getdensa(int *dim = nullptr) const;
  virtual const double *getdensb(int *dim = nullptr) const;
  virtual const double *getfocka(int *dim = nullptr) const;
  virtual const double *getfockb(int *dim = nullptr) const;

  virtual void setcoefa(double *arr);
  virtual void setcoefb(double *arr);
  virtual void setdensa(double *arr);
  virtual void setdensb(double *arr);
  virtual void setfocka(double *arr);
  virtual void setfockb(double *arr);
  
  virtual double getenergy() const;
  virtual const double *getenergies(int *dim) const;
  virtual int getnorbs() const override;

  virtual void setenergies(double *arr);
};

class SCF : public Energy {
public:
  SCF() : Energy() {;}
  SCF(OptionList &opts) : Energy(opts) {;}

  virtual ~SCF() = default;
};

}


#endif
