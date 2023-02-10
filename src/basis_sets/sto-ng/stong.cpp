/*
 * stong.cpp
 *
 *  Created on: Feb 7, 2023
 *      Author: connor
 */

#ifndef __STONG_CPP__
#define __STONG_CPP__

#include "sto3g.hpp"
#include "../electron_configs.hpp"
#include <array>

template<int n> compchem::STO_nGComputer<n>::STO_nGComputer() {
	;
}


template<int n> compchem::STO_nGComputer<n> &compchem::STO_nGComputer<n>::getSingleton() {
	if(compchem::STO_nGComputer<n>::singleton == nullptr) {
		compchem::STO_nGComputer<n>::singleton = new compchem::STO_nGComputer<n>();
	}
	return *compchem::STO_nGComputer<n>::singleton;
}

template<int n> int compchem::STO_nGComputer<n>::getGauss() const {
	return n;
}

static double slater_rule(int n, int l, const compchem::GSConfig &conf) {
	static std::array<int, 19> lvals = std::array<int, 19>({0, 0, 1, 0, 1, 0, 2, 1, 0, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1});
	static std::array<int, 19> nvals = std::array<int, 19>({1, 2, 2, 3, 3, 4, 3, 4, 5, 4, 5, 6, 4, 5, 6, 7, 5, 6, 7});
	if(n == 1) {
		return 0.3 * (conf[0] - 1);
	} else {
		int curr = 0;
		double out = 0;

		for(int i = 0; i < 19; i++) {
			if(nvals[curr] == n && l == 2 && lvals[curr] < l) {
				out += conf[curr];
			} else if(nvals[curr] == n && (l == 0 || l == 1) && (lvals[curr] == 0 || lvals[curr] == 1)) {
				if(lvals[curr] == l) {
					out += 0.35 * (conf[curr] - 1);
				} else {
					out += 0.35 * conf[curr];
				}
			} else if(nvals[curr] == n && lvals[curr] == l) {
				out += 0.35 * (conf[curr] - 1);
			} else if(nvals[curr] == n - 1) {
				if(l == 0 || l == 1) {
					out += 0.85 * conf[curr];
				} else {
					out += conf[curr];
				}
			} else if(nvals[curr] <= n - 2) {
				out += conf[curr];
			}

			curr++;
		}
		return out;
	}
}

template<int n> compchem::GaussianSum *compchem::STO_nGComputer<n>::compute(int Z, int *nouts) {
	static std::array<int, 19> lvals = std::array<int, 19>({0, 0, 1, 0, 1, 0, 2, 1, 0, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1});
	static std::array<int, 19> nvals = std::array<int, 19>({1, 2, 2, 3, 3, 4, 3, 4, 5, 4, 5, 6, 4, 5, 6, 7, 5, 6, 7});

	int orbs = 0;

	for(int i = 0; i < compchem::gsconfigs[Z - 1].getShells(); i++) {
		orbs += 2 * lvals[i] + 1;
	}

	*nouts = orbs;

	compchem::GaussianSum out[] = new GaussianSum[orbs];

	for(int i = 0; i < compchem::gsconfigs[Z - 1].getShells(); i++)  {
		double zeff = Z - slater_rule(nvals[i], lvals[i], compchem::gsconfigs[Z - 1]);
		for(int m = -lvals[i]; m <= lvals[i]; m++) {

		}
	}

}

#endif


