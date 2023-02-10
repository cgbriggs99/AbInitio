/*
 * compute.cpp
 *
 *  Created on: Feb 7, 2023
 *      Author: connor
 */

#include "compute.hpp"
#include <cmath>
#include <cstring>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <array>

using namespace compchem;

int shells(int Z) {
	// TODO: find a closed form perhaps.
	int shell = 0;

	while(Z > 0) {
		shell++;
		Z -= 2 * (shell / 2 + 1) * (shell / 2 + 1);
	}
	return shell;
}

ShellFactory::ShellFactory() {
	;
}

ShellFactory &ShellFactory::getSingleton() {
	if(ShellFactory::singleton == NULL) {
		ShellFactory::singleton = new ShellFactory();
	}
	return *ShellFactory::singleton;
}

GaussianSum &ShellFactory::makeShell(int l, int m, double alpha) {
	GaussianSum *out = new GaussianSum(0);

	std::vector<std::array<int, 3> > *expos1 = new std::vector<std::array<int, 3> >(0);
	std::vector<double> *coefs1 = new std::vector<double>(0);

	std::vector<std::array<int, 3> > *expos2 = new std::vector<std::array<int, 3> >(0);
	std::vector<double> *coefs2 = new std::vector<double>(0);

	if(m < 0) {
		// First, compute Pi(z).
		double fact = exp(lgamma(l + m + 1) / 2 - lgamma(l - m + 1) / 2 + lgamma(2 * l + 1) -
				lgamma(l - m + 1)) / pow(2, l);
		for(int i = 0; i <= (l + m) / 2; i++) {
			// Factor out the r^2k z^(l-2k+m) term.
			for(int j = 0; j <= i; j++) {
				for(int k = 0; k <= i - j; k++) {
					expos1->push_back(*new std::array<int, 3>({2 * i, 2 * j, l - m - 2 * j - 2 * k}));
					coefs1->push_back(fact * exp(lgamma(i + 1) - lgamma(j + 1) - lgamma(k + 1) - lgamma(i - j - k + 1)));
				}
			}
			// Update the coefficient.
			fact *= -(l - i) * (l + m - 2 * i) * (l + m - 2 * i - 1);
			fact /= (i + 1) * (2 * l - 2 * i) * (2 * l - 2 * i - 1);
		}


		// Then, compute B and multiply.
		fact = 1;
		for(int i = 0; i <= -m; i++) {
			if((-m - i) % 4 == 1) {
				for(int j = 0; j < expos1->size(); i++) {
					std::array<int, 3> &expo = expos1->at(j);
					double coef = coefs1->at(j);

					expos2->push_back(*new std::array<int, 3>({expo[0] + i, expo[1] - m - i, expo[2]}));
					coefs2->push_back(coef * fact);
				}
			} else if((-m - i) % 4 == 3) {
				for(int j = 0; j < expos1->size(); i++) {
					std::array<int, 3> &expo = expos1->at(j);
					double coef = coefs1->at(j);

					expos2->push_back(*new std::array<int, 3>({expo[0] + i, expo[1] - m - i, expo[2]}));
					coefs2->push_back(-coef * fact);
				}
			}
			fact *= (-m - i);
			fact /= (i + 1);
		}

		// Finally, form the output.
		for(int i = 0; i < expos2->size(); i++) {
			int a = expos2->at(i)[0], b = expos2->at(i)[1], c = expos2->at(i)[2];
			double coef = coefs2->at(i);
			*out += *new GaussianTerm(coef * pow(2 * alpha / M_PI, 0.75) * sqrt(pow(8 * alpha, l) *
					exp(lgamma(a + 1) + lgamma(b + 1) + lgamma(c + 1)
							- lgamma(2 * a + 1) - lgamma(2 * b + 1) - lgamma(2 * c + 1))),
							std::array<double, 3>({0, 0, 0}), std::array<double, 3>({alpha, alpha, alpha}),
							expos2->at(i));
		}
		delete expos1, expos2, coefs1, coefs2;
		return *out;

	} else if(m > 0) {
		// First, compute Pi(z).
		double fact = exp(lgamma(l - m + 1) / 2 - lgamma(l + m + 1) / 2 + lgamma(2 * l + 1) -
				lgamma(l + m + 1)) / pow(2, l);
		for(int i = 0; i <= (l - m) / 2; i++) {
			// Factor out the r^2k z^(l-2k+m) term.
			for(int j = 0; j <= i; j++) {
				for(int k = 0; k <= i - j; k++) {
					expos1->push_back(*new std::array<int, 3>({2 * i, 2 * j, l + m - 2 * j - 2 * k}));
					coefs1->push_back(fact * exp(lgamma(i + 1) - lgamma(j + 1) - lgamma(k + 1) - lgamma(i - j - k + 1)));
				}
			}
			// Update the coefficient.
			fact *= -(l - i) * (l - m - 2 * i) * (l - m - 2 * i - 1);
			fact /= (i + 1) * (2 * l - 2 * i) * (2 * l - 2 * i - 1);
		}


		// Then, compute B and multiply.
		fact = 1;
		for(int i = 0; i <= m; i++) {
			if((m - i) % 4 == 1) {
				for(int j = 0; j < expos1->size(); i++) {
					std::array<int, 3> &expo = expos1->at(j);
					double coef = coefs1->at(j);

					expos2->push_back(*new std::array<int, 3>({expo[0] + i, expo[1] + m - i, expo[2]}));
					coefs2->push_back(coef * fact);
				}
			} else if((m - i) % 4 == 3) {
				for(int j = 0; j < expos1->size(); i++) {
					std::array<int, 3> &expo = expos1->at(j);
					double coef = coefs1->at(j);

					expos2->push_back(*new std::array<int, 3>({expo[0] + i, expo[1] + m - i, expo[2]}));
					coefs2->push_back(-coef * fact);
				}
			}
			fact *= (m - i);
			fact /= (i + 1);
		}

		// Finally, form the output.
		for(int i = 0; i < expos2->size(); i++) {
			int a = expos2->at(i)[0], b = expos2->at(i)[1], c = expos2->at(i)[2];
			double coef = coefs2->at(i);
			*out += *new GaussianTerm(coef * pow(2 * alpha / M_PI, 0.75) * sqrt(pow(8 * alpha, l) *
					exp(lgamma(a + 1) + lgamma(b + 1) + lgamma(c + 1)
							- lgamma(2 * a + 1) - lgamma(2 * b + 1) - lgamma(2 * c + 1))),
							std::array<double, 3>({0, 0, 0}), std::array<double, 3>({alpha, alpha, alpha}),
							expos2->at(i));
		}
		delete expos1, expos2, coefs1, coefs2;
		return *out;
	} else {
		// First, compute Pi(z).
		double fact = exp(lgamma(2 * l + 1) - lgamma(l + 1)) / pow(2, l);
		for(int i = 0; i <= l / 2; i++) {
			// Factor out the r^2k z^(l-2k+m) term.
			for(int j = 0; j <= i; j++) {
				for(int k = 0; k <= i - j; k++) {
					*out += *new GaussianTerm(fact *
							exp(lgamma(i + 1) - lgamma(j + 1) - lgamma(k + 1) - lgamma(i - j - k + 1)),
							std::array<double,3>({0, 0, 0}), std::array<double, 3>({alpha, alpha, alpha}),
							std::array<int, 3>({2 * i, 2 * j, l - m - 2 * j - 2 * k}));
				}
			}
			// Update the coefficient.
			fact *= -(l - i) * (l + m - 2 * i) * (l + m - 2 * i - 1);
			fact /= (i + 1) * (2 * l - 2 * i) * (2 * l - 2 * i - 1);
		}

		// Finally, form the output.
		for(int i = 0; i < expos2->size(); i++) {
			int a = expos2->at(i)[0], b = expos2->at(i)[1], c = expos2->at(i)[2];
			double coef = coefs2->at(i);
			*out += *new GaussianTerm(coef * pow(2 * alpha / M_PI, 0.75) * sqrt(pow(8 * alpha, l) *
					exp(lgamma(a + 1) + lgamma(b + 1) + lgamma(c + 1)
							- lgamma(2 * a + 1) - lgamma(2 * b + 1) - lgamma(2 * c + 1))),
							std::array<double, 3>({0, 0, 0}), std::array<double, 3>({alpha, alpha, alpha}),
							expos2->at(i));
		}
		delete expos1, expos2, coefs1, coefs2;
		return *out;
	}
}
