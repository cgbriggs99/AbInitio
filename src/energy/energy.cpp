#include "energy.hpp"


using namespace compchem;

OptionList &Energy::getoptions() {
  return this->options;
}

const OptionList &Energy::getoptions() const {
  return this->options;
}

