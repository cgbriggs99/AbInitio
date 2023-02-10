/*
 * integrals.h
 *
 *  Created on: Feb 6, 2023
 *      Author: connor
 */

#ifndef INTEGRALS_HPP_
#define INTEGRALS_HPP_

#include "gaussian.hpp"

namespace compchem {

double gaussianintegral(const GaussianSum &overlap);

double oeintegral(const GaussianSum &overlap, const std::array<double, 3> center);

double teintegral(const GaussianSum &electron1, const GaussianSum &electron2);

};

#endif /* INTEGRALS_HPP_ */
