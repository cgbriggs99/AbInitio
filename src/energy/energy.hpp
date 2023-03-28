#ifndef ENERGY_HPP
#define ENERGY_HPP

#include "../wavefunction/wavefunction.hpp"
#include "../opts/options.hpp"
#include "../util/molecule.hpp"
#include "../integrals/integrals.hpp"

namespace compchem {

class Energy {
protected:
  OptionList &options;

public:
  Energy() : options(GlobalOptions::getsingleton()) {;}
  Energy(OptionList &opts) : options(opts) {;}

  virtual ~Energy() = default;

  virtual OptionList &getoptions();
  virtual const OptionList &getoptions() const;

  virtual double energy(const Molecule *molecule,
			const Wavefunction *wfn_in) const = 0;
};

}

#endif
