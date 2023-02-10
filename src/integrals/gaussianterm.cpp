/*
 * gaussian.cpp
 *
 *  Created on: Nov 8, 2022
 *      Author: connor
 */

#include "gaussian.hpp"
#include <array>
#include <cmath>
#include <stdexcept>

using namespace compchem;

GaussianTerm::GaussianTerm(double coef, const std::array<double, 3> &center, const std::array<double, 3> &mults,
		const std::array<int, 3> &expos) {
	this->coef = coef;
	this->center = center;
	this->mults = mults;
	this->expos = expos;
	this->scalar = (mults[0] == 0 && mults[1] == 0 && mults[2] == 0 &&
			expos[0] == 0 && expos[1] == 0 && expos[2] == 0) || coef == 0;
}

GaussianTerm::GaussianTerm(double cast) {
	this->center.fill(0);
	this->coef = cast;
	this->expos.fill(0);
	this->mults.fill(0);
	this->scalar = true;
}

GaussianTerm::GaussianTerm() {
	this->center.fill(0);
	this->coef = 0;
	this->expos.fill(0);
	this->mults.fill(0);
	this->scalar = true;
}

GaussianTerm::GaussianTerm(const GaussianTerm &copy) {
	this->coef = copy.coef;
	this->scalar = copy.scalar;
	for(int i = 0; i < 3; i++) {
		this->center[i] = copy.center[i];
		this->expos[i] = copy.expos[i];
		this->mults[i] = copy.mults[i];
	}
}

//GaussianTerm::GaussianTerm(GaussianTerm &&move) {
//	this->coef = move.coef;
//	for(int i = 0; i < 3; i++) {
//		this->center[i] = move.center[i];
//		this->expos[i] = move.expos[i];
//		this->mults[i] = move.mults[i];
//	}
//}

double GaussianTerm::getCoef() const {
	return this->coef;
}

const std::array<double, 3> &GaussianTerm::getCenter() const {
	return this->center;
}

const std::array<int, 3> &GaussianTerm::getExpos() const {
	return this->expos;
}

const std::array<double, 3> &GaussianTerm::getMults() const {
	return this->mults;
}

double GaussianTerm::getCenter(int dim) const {
	if(dim < 0 || dim >= 3) {
		throw new std::out_of_range("Dimension not between 0 and 3.");
	}
	return this->center[dim];
}

int GaussianTerm::getExpo(int dim) const {
	if(dim < 0 || dim >= 3) {
		throw new std::out_of_range("Dimension not between 0 and 3.");
	}
	return this->expos[dim];
}

double GaussianTerm::getMult(int dim) const {
	if(dim < 0 || dim >= 3) {
		throw new std::out_of_range("Dimension not between 0 and 3.");
	}
	return this->mults[dim];
}

void GaussianTerm::setCoef(double coef) {
	this->coef = coef;
}

void GaussianTerm::setCenter(const std::array<double, 3> &center) {
	for(int i = 0; i < 3; i++) {
		this->center[i] = center[i];
	}
}

void GaussianTerm::setExpos(const std::array<int, 3> &expos) {
	for(int i = 0; i < 3; i++) {
		this->expos[i] = expos[i];
	}
}

void GaussianTerm::setMults(const std::array<double, 3> &mults) {
	for(int i = 0; i < 3; i++) {
		this->mults[i] = mults[i];
	}
}

double GaussianTerm::eval(double x, double y, double z) const {
	return this->coef * pow(x, this->expos[0]) * pow(y, this->expos[1]) * pow(z, this->expos[2]) *
			exp(this->mults[0] * pow(x - this->center[0], 2) + this->mults[1] * pow(y - this->center[1], 2) +
					this->mults[2] * pow(z - this->center[2], 2));
}

double GaussianTerm::eval(std::array<double, 3> coord) const {
	return this->coef * pow(coord[0], this->expos[0]) * pow(coord[1], this->expos[1]) * pow(coord[2], this->expos[2]) *
				exp(this->mults[0] * pow(coord[0] - this->center[0], 2) + this->mults[1] * pow(coord[1] - this->center[1], 2) +
						this->mults[2] * pow(coord[2] - this->center[2], 2));
}

bool GaussianTerm::orderEqual(const GaussianTerm &other) const {
	return (this->expos[0] == other.expos[0] && this->expos[1] == other.expos[1] && this->expos[2] == other.expos[2] &&
			this->center[0] == other.center[0] && this->center[1] == other.center[1] && this->center[2] == other.center[2] &&
			this->mults[0] == other.mults[0] && this->mults[1] == other.mults[1] && this->mults[2] == other.mults[2]) ||
			(this->scalar == other.scalar);
}

bool GaussianTerm::isScalar() const {
	return this->scalar;
}

GaussianTerm &GaussianTerm::operator*=(const GaussianTerm &rh) {
	/*
	 * a (x - b)^2 + c (x - d)^2
	 * ax^2 - 2axb + ab^2 + cx^2 - 2cdx + cd^2
	 * (a + c) x^2 - 2(ab + cd) x + ab^2 + cd^2
	 * y = (ab + cd) / (a + c)
	 * z = (ab^2 + cd^2) / (a + c)
	 * (a + c) (x^2 - 2xy + z)
	 * (a + c) (x^2 - 2xy + y^2 - y^2 + z)
	 * (a + c) (x - y)^2 + (a + c) (z - y^2)
	 */
	double ax = this->getMult(0), ay = this->getMult(1), az = this->getMult(2);
	double bx = this->getCenter(0), by = this->getCenter(1), bz = this->getCenter(2);
	double cx = rh.getMult(0), cy = rh.getMult(1), cz = rh.getMult(2);
	double dx = rh.getCenter(0), dy = rh.getCenter(1), dz = rh.getCenter(2);

	double centx = (ax * bx + cx * dx) / (ax + cx), centy = (ay * by + cy * dy) / (ay + cy),
			centz = (az * bz + cz * dz) / (az + cz);
	double corrx = (ax * bx * bx + cx * dx * dx) / (ax + cx), corry = (ay * by * by + cy * dy * dy) / (ay + cy),
			corrz = (az * bz * bz + cz * dz * dz) / (az + cz);

	this->center[0] = (ax * bx + cx * dx) / (ax + cx);
	this->center[1] = (ay * by + cy * dy) / (ay + cy);
	this->center[2] = (az * bz + cz * dz) / (az + cz);

	this->mults[0] = ax + cx;
	this->mults[1] = ay + cy;
	this->mults[2] = az + cz;

	this->expos[0] += rh.getExpo(0);
	this->expos[1] += rh.getExpo(1);
	this->expos[2] += rh.getExpo(2);

	this->coef *= rh.coef * exp((ax + cx) * (corrx - centx * centx) + (ay + cy) * (corry - centy * centy) +
			(az + cz) * (corrz - centz * centz));

	this->scalar = (this->isScalar() && rh.isScalar()) || (this->mults[0] == 0 && this->mults[1] == 0 && this->mults[2] == 0
			&& this->expos[0] == 0 && this->expos[1] == 0 && this->expos[2] == 0) || (this->coef == 0);

	return *this;
}

GaussianTerm &GaussianTerm::operator*=(double scal) {
	this->coef *= scal;
	return *this;
}

GaussianTerm &GaussianTerm::operator/=(double scal) {
	this->coef /= scal;
	return *this;
}

GaussianTerm &GaussianTerm::operator*(const GaussianTerm &rh) const {
	GaussianTerm *out = new GaussianTerm(*this);

	*out *= rh;
	return *out;
}

GaussianTerm &GaussianTerm::operator*(double rh) const {
	GaussianTerm *out = new GaussianTerm(*this);

	*out *= rh;
	return *out;
}

GaussianTerm &GaussianTerm::operator/(double rh) const {
	GaussianTerm *out = new GaussianTerm(*this);

	*out /= rh;
	return *out;
}

GaussianTerm &GaussianTerm::operator-() const {
	GaussianTerm *out = new GaussianTerm(*this);

	out->coef = -out->coef;
	return *out;
}

GaussianSum &GaussianTerm::operator+(const GaussianTerm &rh) {
	if(this->orderEqual(rh)) {
		return GaussianSum(*this) + rh;
	}
	return *new GaussianSum(2, this, &rh);
}

GaussianSum &GaussianTerm::operator+(const GaussianSum &rh) {
	return rh + *this;
}

GaussianSum &GaussianTerm::operator+(double rh) {
	return GaussianSum(*this) + rh;
}

GaussianSum &GaussianTerm::operator-(const GaussianTerm &rh) {
	if(this->orderEqual(rh)) {
		return GaussianSum(*this) - rh;
	}
	return *new GaussianSum(2, this, -rh);
}

GaussianSum &GaussianTerm::operator-(const GaussianSum &rh) {
	return -rh + *this;
}

GaussianSum &GaussianTerm::operator-(double rh) {
	return *this + (-rh);
}

GaussianSum &GaussianTerm::operator*(const GaussianSum &rh) {
	return rh * *this;
}

