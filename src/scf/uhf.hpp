#ifndef UHF_HPP
#define UHF_HPP

#include "scf.hpp"

namespace compchem {

class UHFWfn : public SCFWfn {
public :
  UHFWfn(int electrons, int orbs, double *S = nullptr, double *T = nullptr,
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
  virtual ~UHFWfn() = default;

  virtual const double *getcoefa(int *dim = nullptr) const override;
  virtual const double *getcoefb(int *dim = nullptr) const override;
  virtual const double *getdensa(int *dim = nullptr) const override;
  virtual const double *getdensb(int *dim = nullptr) const override;
  virtual const double *getfocka(int *dim = nullptr) const override;
  virtual const double *getfockb(int *dim = nullptr) const override;

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

class UHF : public SCF {
public:
  UHF() : SCF() {;}
  UHF(OptionList &opts) : SCF(opts) {;}

  virtual ~UHF() = default;

  virtual double energy(const Molecule *molecule,
			const Wavefunction *wfn_in) const override;
  virtual double energy(const Molecule *molecule,
			const Wavefunction *wfn_in,
		        UHFWfn *wfn_out,
			std::vector<int> occupieda = std::vector<int>{},
			std::vector<int> occupiedb = std::vector<int>{}) const;
};

}


#endif
