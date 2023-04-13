#ifndef RHF_HPP
#define RHF_HPP

#include "scf.hpp"

namespace compchem {

class RHFWfn : public SCFWfn {
public :
  RHFWfn(int electrons, int orbs, double *S = nullptr, double *T = nullptr,
	 double *V = nullptr,
	 TEIArray *tei = nullptr, double *C = nullptr,
	 double *D = nullptr,
	 double *F = nullptr,
	 double *es = nullptr,
	 int multiplicity = 1) : SCFWfn(electrons, orbs, S, T, V, tei,
					C, nullptr,
					D, nullptr,
					F, nullptr,
					es, multiplicity) {;}
  virtual ~RHFWfn() = default;

  virtual const double *getcoefa(int *dim = nullptr) const override;
  virtual const double *getcoefb(int *dim = nullptr) const override;
  virtual const double *getdensa(int *dim = nullptr) const override;
  virtual const double *getdensb(int *dim = nullptr) const override;
  virtual const double *getfocka(int *dim = nullptr) const override;
  virtual const double *getfockb(int *dim = nullptr) const override;

  virtual const double *getcoef(int *dim = nullptr) const;
  virtual const double *getdens(int *dim = nullptr) const;
  virtual const double *getfock(int *dim = nullptr) const;

  virtual void setcoefa(double *arr) override;
  virtual void setcoefb(double *arr) override;
  virtual void setdensa(double *arr) override;
  virtual void setdensb(double *arr) override;
  virtual void setfocka(double *arr) override;
  virtual void setfockb(double *arr) override;

  virtual void setcoef(double *arr);
  virtual void setdens(double *arr);
  virtual void setfock(double *arr);
};

class RHF : public SCF {
public:
  RHF() : SCF() {;}
  RHF(OptionList &opts) : SCF(opts) {;}

  virtual ~RHF() = default;

  virtual double energy(const Molecule *molecule,
			const Wavefunction *wfn_in) const override;
  virtual double energy(const Molecule *molecule,
			const Wavefunction *wfn_in,
		        RHFWfn *wfn_out) const;
};

}


#endif
