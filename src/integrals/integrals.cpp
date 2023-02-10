/*
 * integrals.cpp
 *
 *  Created on: Feb 6, 2023
 *      Author: connor
 */


#include <cmath>
#include "integrals.hpp"

using namespace compchem;

double gaussianintegral(const GaussianSum &overlap) {

	double sum = 0;

	for(int i = 0; i < overlap.getNTerms(); i++) {
		if(overlap[i].getCoef() == 0) {
			continue;
		}
		// These terms are odd and therefore have an integral of zero.
		if(overlap[i].getExpo(0) % 2 == 1 || overlap[i].getExpo(1) % 2 == 1 || overlap[i].getExpo(2) % 2 == 1) {
			continue;
		}
		sum += overlap[i].getCoef() / (pow(overlap[i].getMult(0), overlap[i].getExpo(0) / 2) *
				pow(overlap[i].getMult(1), overlap[i].getExpo(1) / 2) *
				pow(overlap[i].getMult(2), overlap[i].getExpo(2) / 2)) *
				tgamma((overlap[i].getExpo(0) + 1.0) / 2) *
				tgamma((overlap[i].getExpo(1) + 1.0) / 2) *
				tgamma((overlap[i].getExpo(2) + 1.0) / 2);
	}

	return sum;

}

static int num_points = 30;
static double hermite_roots[] = {
		-6.86334529353356,
		6.86334529353359,
		-6.13827922010617,
		6.13827922010600,
		-5.53314715160243,
		5.53314715160413,
		-4.98891896855491,
		4.98891896854903,
		-4.48305535711061,
		4.48305535712055,
		-4.00390860385509,
		4.00390860384527,
		-3.54444387316212,
		3.54444387316849,
		-3.09997052957755,
		3.09997052957458,
		-2.66713212454226,
		2.66713212454325,
		-2.24339146775844,
		2.24339146775825,
		-1.82674114360458,
		1.82674114360458,
		-1.41552780019803,
		1.41552780019804,
		-1.00833827104674,
		1.00833827104674,
		-0.603921058625551,
		0.603921058625552,
		-0.201128576548871,
		0.201128576548872
};

static double hermite_weights[] = {
		2.90825469983025e-21,
		2.90825469983435e-21,
		2.81033360394034e-17,
		2.81033360390304e-17,
		2.87860707840174e-14,
		2.87860707812132e-14,
		8.10618630078503e-12,
		8.10618630891581e-12,
		9.17858041870959e-10,
		9.17858041191871e-10,
		5.10852245172912e-08,
		5.10852245423969e-08,
		1.57909488757581e-06,
		1.57909488706843e-06,
		2.93872522900291e-05,
		2.93872522912957e-05,
		0.000348310124309906,
		0.000348310124278314,
		0.00273792247319975,
		0.00273792247316205,
		0.0147038297047570,
		0.0147038297047247,
		0.0551441768702727,
		0.0551441768702696,
		0.146735847540883,
		0.146735847540888,
		0.280130930839211,
		0.280130930839212,
		0.386394889541814,
		0.386394889541814
};

double oeintegral(const GaussianSum &overlap, const std::array<double, 3> center) {
	double sum = 0;

	for(int i = 0; i < overlap.getNTerms(); i++) {
		for(int x = 0; x < num_points; x++) {
			for(int y = 0; y < num_points; y++) {
				for(int z = 0; z < num_points; z++) {
					sum += hermite_weights[x] * hermite_weights[y] * hermite_weights[z] * overlap[i].getCoef() *
							pow(hermite_roots[x] / sqrt(overlap[i].getMult(0)), overlap[i].getExpo(0)) *
							pow(hermite_roots[y] / sqrt(overlap[i].getMult(1)), overlap[i].getExpo(1)) *
							pow(hermite_roots[z] / sqrt(overlap[i].getMult(2)), overlap[i].getExpo(2)) /
							(sqrt(overlap[i].getMult(0) * overlap[i].getMult(1) * overlap[i].getMult(2)) *
									sqrt((center[0] - hermite_roots[x]) * (center[0] - hermite_roots[x]) +
											(center[1] - hermite_roots[y]) * (center[1] - hermite_roots[y]) +
											(center[2] - hermite_roots[z]) * (center[2] - hermite_roots[z])));
				}
			}
		}
	}
	return sum;
}

double teintegral(const GaussianSum &electron1, const GaussianSum &electron2) {
	double sum = 0;

	for(int i = 0; i < electron1.getNTerms(); i++) {
		for(int j = 0; j < electron2.getNTerms(); j++) {
			for(int x1 = 0; x1 < num_points; x1++) {
				for(int y1 = 0; y1 < num_points; y1++) {
					for(int z1 = 0; z1 < num_points; z1++) {
						for(int x2 = 0; x2 < num_points; x2++) {
							for(int y2 = 0; y2 < num_points; y2++) {
								for(int z2 = 0; z2 < num_points; z2++) {
									sum += hermite_weights[x1] * hermite_weights[y1] * hermite_weights[z1] *
											hermite_weights[x2] * hermite_weights[y2] * hermite_weights[z2] *
											electron1[i].getCoef() * electron2[j].getCoef() *
											pow(hermite_roots[x1] / sqrt(electron1[i].getMult(0)), electron1[i].getExpo(0)) *
											pow(hermite_roots[y1] / sqrt(electron1[i].getMult(1)), electron1[i].getExpo(1)) *
											pow(hermite_roots[z1] / sqrt(electron1[i].getMult(2)), electron1[i].getExpo(2)) *
											pow(hermite_roots[x2] / sqrt(electron2[j].getMult(0)), electron2[j].getExpo(0)) *
											pow(hermite_roots[y2] / sqrt(electron2[j].getMult(1)), electron2[j].getExpo(1)) *
											pow(hermite_roots[z2] / sqrt(electron2[j].getMult(2)), electron2[j].getExpo(2))/
											(sqrt(electron1[i].getMult(0) * electron1[i].getMult(1) * electron1[i].getMult(2) *
													electron2[j].getMult(0) * electron2[j].getMult(1) * electron2[j].getMult(2)) *
													sqrt((hermite_roots[x1] - hermite_roots[x2]) *
															(hermite_roots[x1] - hermite_roots[x2]) +
															(hermite_roots[y1] - hermite_roots[y2]) *
															(hermite_roots[y1] - hermite_roots[y2]) +
															(hermite_roots[z1] - hermite_roots[z2]) *
															(hermite_roots[z1] - hermite_roots[z2])));
								}
							}
						}
					}
				}
			}
		}
	}
	return sum;
}

