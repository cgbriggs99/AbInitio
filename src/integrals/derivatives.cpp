/*
 * derivatives.cpp
 *
 *  Created on: Nov 8, 2022
 *      Author: connor
 */



#include "derivatives.hpp"
#include "gaussian.hpp"

using namespace compchem;

GaussianSum &laplacian(const GaussianTerm &term) {
	/*
	 * x^n e^(a(x - b)^2)
	 * nx^(n-1) e^(a(x - b)^2) + 2a(x - b)x^n e^(a(x - b)^2)
	 * n(n-1)x^(n-2) e^(a(x - b)^2) + 2na(x-b)x^(n-1) e^(a(x - b)^2) + 2ax^n e^(a(x - b)^2) +
	 *  2na(x - b)x^(n-1) e^(a(x - b)^2) + 4a^2(x - b)^2x^n e^(a(x - b)^2)
	 * (n(n-1)x^(n-2) + 4na(x-b)x^(n-1) + 2ax^n + 4a^2(x-b)^2x^n)e^(a(x-b)^2)
	 * (n(n-1)x^(n-2) - 4nabx^(n-1) + (4na + 2a + 4a^2b^2)x^n - 8a^2bx^(n+1) + 4a^2x^(n+2))e^(a(x-b)^2)
	 */

	std::array<int, 3> px1, px2, px3, px4, py1, py2, py3, py4, pz1, pz2, pz3, pz4;
	px1[0] = term.getExpo(0) - 2, px1[1] = term.getExpo(1), px1[2] = term.getExpo(2);
	px2[0] = term.getExpo(0) - 1, px2[1] = term.getExpo(1), px2[2] = term.getExpo(2);
	px3[0] = term.getExpo(0) + 1, px3[1] = term.getExpo(1), px3[2] = term.getExpo(2);
	px4[0] = term.getExpo(0) + 2, px4[1] = term.getExpo(1), px4[2] = term.getExpo(2);
	py1[0] = term.getExpo(0), py1[1] = term.getExpo(1) - 2, py1[2] = term.getExpo(2);
	py2[0] = term.getExpo(0), py2[1] = term.getExpo(1) - 1, py2[2] = term.getExpo(2);
	py3[0] = term.getExpo(0), py3[1] = term.getExpo(1) + 1, py3[2] = term.getExpo(2);
	py4[0] = term.getExpo(0), py4[1] = term.getExpo(1) + 2, py4[2] = term.getExpo(2);
	pz1[0] = term.getExpo(0), pz1[1] = term.getExpo(1), pz1[2] = term.getExpo(2) - 2;
	pz2[0] = term.getExpo(0), pz2[1] = term.getExpo(1), pz2[2] = term.getExpo(2) - 1;
	pz3[0] = term.getExpo(0), pz3[1] = term.getExpo(1), pz3[2] = term.getExpo(2) + 1;
	pz4[0] = term.getExpo(0), pz4[1] = term.getExpo(1), pz4[2] = term.getExpo(2) + 2;

	return term.getCoef() * (GaussianTerm(term.getExpo(0) * (term.getExpo(0) - 1), term.getCenter(), term.getMults(), px1) -
			GaussianTerm(4 * term.getExpo(0) * term.getMult(0) * term.getCenter(0),
					term.getCenter(), term.getMults(), px2) +
			GaussianTerm(4 * term.getExpo(0) * term.getMult(0) + 2 * term.getMult(0) +
					4 * term.getMult(0) * term.getMult(0) * term.getCenter(0) * term.getCenter(0),
					term.getCenter(), term.getMults(), term.getExpos()) -
			GaussianTerm(8 * term.getMult(0) * term.getMult(0) * term.getCenter(0), term.getCenter(), term.getMults(), px3) +
			GaussianTerm(4 * term.getMult(0) * term.getMult(0), term.getCenter(), term.getMults(), px4) +
			GaussianTerm(term.getExpo(1) * (term.getExpo(1) - 1), term.getCenter(), term.getMults(), py1) -
			GaussianTerm(4 * term.getExpo(1) * term.getMult(1) * term.getCenter(1),
					term.getCenter(), term.getMults(), py2) +
			GaussianTerm(4 * term.getExpo(1) * term.getMult(1) + 2 * term.getMult(1) +
					4 * term.getMult(1) * term.getMult(1) * term.getCenter(1) * term.getCenter(1),
					term.getCenter(), term.getMults(), term.getExpos()) -
			GaussianTerm(8 * term.getMult(1) * term.getMult(1) * term.getCenter(1), term.getCenter(), term.getMults(), py3) +
			GaussianTerm(4 * term.getMult(1) * term.getMult(1), term.getCenter(), term.getMults(), py4) +
			GaussianTerm(term.getExpo(2) * (term.getExpo(2) - 1), term.getCenter(), term.getMults(), pz1) -
			GaussianTerm(4 * term.getExpo(2) * term.getMult(2) * term.getCenter(2),
					term.getCenter(), term.getMults(), pz2) +
			GaussianTerm(4 * term.getExpo(2) * term.getMult(2) + 2 * term.getMult(2) +
					4 * term.getMult(2) * term.getMult(2) * term.getCenter(2) * term.getCenter(2),
					term.getCenter(), term.getMults(), term.getExpos()) -
			GaussianTerm(8 * term.getMult(2) * term.getMult(2) * term.getCenter(2), term.getCenter(), term.getMults(), pz3) +
			GaussianTerm(4 * term.getMult(2) * term.getMult(2), term.getCenter(), term.getMults(), pz4));
}

GaussianSum &laplacian(const GaussianSum &sum) {
	GaussianSum *out = new GaussianSum(0);

	for(int i = 0; i < sum.getNTerms(); i++) {
		*out += compchem::laplacian(sum[i]);
	}
	return *out;
}
