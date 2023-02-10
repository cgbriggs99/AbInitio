/*
 * sto3g.hpp
 *
 *  Created on: Feb 7, 2023
 *      Author: connor
 */

#ifndef STO3G_HPP_
#define STO3G_HPP_

#include "../compute.hpp"

namespace compchem {
template<int n>
class STO_nGComputer : OrbitalComputer {
private :
	static STO_nGComputer *singleton = NULL;

	STO_nGComputer();

public:
	~STO_nGComputer() = default;

	static STO_nGComputer<n> &getSingleton();

	int getGauss() const;

	GaussianSum *compute(int Z, int *norbs) override;

};

// For simplicity.
typedef class STO_3GComputer STO_nGComputer<3>;

}

#include "stong.cpp"

#endif /* STO3G_HPP_ */
