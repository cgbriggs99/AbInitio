#include "options.hpp"
#include "default_options.hpp"

using namespace compchem;
using namespace std;

void DefaultOptionsFactory::initializeoptions(OptionList &opts) {

  opts.setintoption("threads", 4);
  opts.setbooloption("analytic boys", true);
  opts.setintoption("boys points", 32);
  opts.setintoption("max scf cycles", 100);
  opts.setdoubleoption("scf rms convergence", 1e-6);
  opts.setdoubleoption("scf energy convergence", 1e-7);

}
