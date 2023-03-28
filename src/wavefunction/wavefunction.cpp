#include "wavefunction.hpp"

using namespace compchem;

int Wavefunction::getelectrons() const {
  return this->electrons;
}

int Wavefunction::getmultiplicity() const {
  return this->multiplicity;
}
