/*
 * derivatives.hpp
 *
 *  Created on: Nov 8, 2022
 *      Author: connor
 */

#ifndef DERIVATIVES_HPP_
#define DERIVATIVES_HPP_

#include "gaussian.hpp"

namespace compchem {

GaussianSum &laplacian(const GaussianTerm &term);
GaussianSum &laplacian(const GaussianSum &sum);

}

#endif /* DERIVATIVES_HPP_ */
