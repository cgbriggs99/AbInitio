/*
 * ints.hpp
 *
 *  Created on: Nov 8, 2022
 *      Author: connor
 */

#ifndef GAUSSIAN_HPP_
#define GAUSSIAN_HPP_

#include <array>
#include <vector>
#include <initializer_list>

namespace compchem {

class Orbital {
public :
	virtual ~Orbital() = 0;

	virtual double eval(double x, double y, double z) const;
};

class GaussianOrbital : Orbital {
private :
	double alpha;
	int l, m;

public :
	GaussianOrbital(double alpha, int l, int ml) :
		alpha(alpha), l(l), m(ml) {
		;
	}

	~GaussianOrbital() = default;

	int getL() const;
	int getM() const;
	double getAlpha() const;

	double eval(double x, double y, double z) const override;
};

class LCGaussianOrbital : Orbital {
private :
	int l, m, norbs;
	double alphas[];
	double coefs[];
public :
	LCGaussianOrbital(int l, int ml, int norbs, double alphas[], double coefs[]);
	LCGaussianOrbital(int l, int ml, int norbs, std::vector<double> alphas, std::vector<double> coefs);
	LCGaussianOrbital(int l, int ml, int norbs, std::initializer_list<double> alphas, std::initializer_list<double> coefs);

	~LCGaussianOrbital();

	double eval(double x, double y, double z) const override;
};

class SlaterOrbital : Orbital {
private :
	double zeff;
	int n, l, m;

public :
	SlaterOrbital(double zeff, int n, int l, int ml) :
		zeff(zeff), n(n), l(l), m(ml) {
		;
	}

	~SlaterOrbital() = default;

	int getN() const;
	int getL() const;
	int getM() const;
	double getZEff() const;

	double operator()(double x, double y, double z) const;

	double eval(double x, double y, double z) const override;
};

}

#endif /* GAUSSIAN_HPP_ */
