/*
 * gaussiansum.cpp
 *
 *  Created on: Nov 8, 2022
 *      Author: connor
 */


#include "gaussian.hpp"
#include <array>
#include <vector>
#include <stdarg.h>
#include <stdexcept>
#include <cmath>

using namespace compchem;

GaussianSum::GaussianSum(int nterms, const GaussianTerm *t1, ...) {
	va_list args;
	va_start(args, t1);

	this->nterms = nterms;
	this->terms = std::vector<GaussianTerm *>(nterms);

	for(int i = 0; i < nterms; i++) {
		this->terms[i] = new GaussianTerm(*va_arg(args, const GaussianTerm *));
	}
	va_end(args);
}

GaussianSum::GaussianSum(const GaussianSum &copy) {
	this->nterms = copy.getNTerms();
	this->terms = std::vector<GaussianTerm *>(this->nterms);

	for(int i = 0; i < this->nterms; i++) {
		this->terms[i] = new GaussianTerm(copy[i]);
	}
}

//
//GaussianSum::GaussianSum(GaussianSum &&move) {
//	this->nterms = move.getNTerms();
//	this->terms(this->nterms);
//
//	for(int i = 0; i < this->nterms; i++) {
//		this->terms[i] = GaussianTerm(move[i]);
//	}
//}

GaussianSum::~GaussianSum() {
	for(int i = 0; i < this->nterms; i++) {
		delete this->terms[i];
	}
}

GaussianSum::GaussianSum(double cast) {
	this->nterms = 1;
	this->terms = std::vector<GaussianTerm*>(1);
	this->terms[0] = new GaussianTerm(cast);
}

GaussianSum::GaussianSum() {
	this->nterms = 1;
	this->terms = std::vector<GaussianTerm*>(1);
	this->terms[0] = new GaussianTerm(0);
}

GaussianSum::GaussianSum(const GaussianTerm &cast) {
	this->nterms = 1;
	this->terms = std::vector<GaussianTerm*>(1);
	this->terms[0] = new GaussianTerm(cast);
}

int GaussianSum::getNTerms() const {
	return this->nterms;
}

GaussianTerm &GaussianSum::operator[](int index) {
	if(index < 0 || index > this->getNTerms()) {
		throw new std::out_of_range("Index out of range in GaussianSum!");
	}
	return *this->terms.at(index);
}

const GaussianTerm &GaussianSum::operator[](int index) const {
	if(index < 0 || index > this->getNTerms()) {
		throw new std::out_of_range("Index out of range in GaussianSum!");
	}
	return *this->terms.at(index);
}

GaussianSum &GaussianSum::operator+=(const GaussianTerm &rh) {
	if(rh.getCoef() == 0) {
		return *this;
	}
	for(int i = 0; i < this->getNTerms(); i++) {
		if(this->terms[i]->orderEqual(rh)) {
			this->terms[i]->setCoef(this->terms[i]->getCoef() + rh.getCoef());
			if(this->terms[i]->getCoef() == 0) {
				this->terms.at(i) = std::move(this->terms.back());
				this->terms.pop_back();
				this->nterms--;
			}
			return *this;
		}
	}

	this->nterms++;
	GaussianTerm copy = GaussianTerm(rh);
	this->terms.push_back(&copy);
	return *this;
}

GaussianSum &GaussianSum::operator+=(const GaussianSum &rh) {
	for(int i = 0; i < rh.getNTerms(); i++) {
		for(int j = 0; j < this->getNTerms(); j++) {
			if(this->terms[j]->orderEqual(rh[i])) {
				this->terms[j]->setCoef(this->terms[j]->getCoef() + rh[i].getCoef());
				if(this->terms[j]->getCoef() == 0) {
					this->terms.at(j) = std::move(this->terms.back());
					this->terms.pop_back();
					this->nterms--;
				}
				break;
			}
		}
		this->nterms++;
		this->terms.push_back(new GaussianTerm(rh[i]));
	}
	return *this;
}

GaussianSum &GaussianSum::operator+=(double rh) {
	if(rh == 0) {
		return *this;
	}
	for(int i = 0; i < this->getNTerms(); i++) {
		if(this->terms[i]->isScalar()) {
			this->terms[i]->setCoef(this->terms[i]->getCoef() + rh);
			if(this->terms[i]->getCoef() == 0) {
				this->terms.at(i) = std::move(this->terms.back());
				this->terms.pop_back();
				this->nterms--;
			}
			return *this;
		}
	}

	this->nterms++;
	this->terms.push_back(new GaussianTerm(rh));
	return *this;
}

GaussianSum &GaussianSum::operator-=(const GaussianTerm &rh) {
	if(rh.getCoef() == 0) {
		return *this;
	}
	for(int i = 0; i < this->getNTerms(); i++) {
		if(this->terms[i]->orderEqual(rh)) {
			this->terms[i]->setCoef(this->terms[i]->getCoef() - rh.getCoef());
			if(this->terms[i]->getCoef() == 0) {
				this->terms.at(i) = std::move(this->terms.back());
				this->terms.pop_back();
				this->nterms--;
			}
			return *this;
		}
	}

	this->nterms++;
	this->terms.push_back(new GaussianTerm(rh));
	return *this;
}

GaussianSum &GaussianSum::operator-=(const GaussianSum &rh) {
	for(int i = 0; i < rh.getNTerms(); i++) {
		for(int j = 0; j < this->getNTerms(); j++) {
			if(this->terms[j]->orderEqual(rh[i])) {
				this->terms[j]->setCoef(this->terms[j]->getCoef() - rh[i].getCoef());
				if(this->terms[j]->getCoef() == 0) {
					this->terms.at(j) = std::move(this->terms.back());
					this->terms.pop_back();
					this->nterms--;
				}
				break;
			}
		}
		this->nterms++;
		this->terms.push_back(new GaussianTerm(rh[i]));
	}
	return *this;
}

GaussianSum &GaussianSum::operator-=(double rh) {
	if(rh == 0) {
		return *this;
	}
	for(int i = 0; i < this->getNTerms(); i++) {
		if(this->terms[i]->isScalar()) {
			this->terms[i]->setCoef(this->terms[i]->getCoef() - rh);
			if(this->terms[i]->getCoef() == 0) {
				this->terms.at(i) = std::move(this->terms.back());
				this->terms.pop_back();
				this->nterms--;
			}
			return *this;
		}
	}

	this->nterms++;
	this->terms.push_back(new GaussianTerm(rh));
	return *this;
}

GaussianSum &GaussianSum::operator*=(const GaussianSum &rh) {
	GaussianSum out = *this * rh;
	*this = std::move(out);
	return *this;
}

GaussianSum &GaussianSum::operator*=(const GaussianTerm &rh) {
	if(rh.getCoef() == 0) {
		this->nterms = 0;
		this->terms.clear();
	}
	for(int i = 0; i < this->getNTerms(); i++) {
		*this->terms[i] *= rh;
		if(this->terms[i]->getCoef() == 0) {
			this->terms.at(i) = std::move(this->terms.back());
			this->terms.pop_back();
			this->nterms--;
			i--;
			continue;
		}
	}
	return *this;
}

GaussianSum &GaussianSum::operator*=(double rh) {
	if(rh == 0) {
		this->nterms = 0;
		this->terms.clear();
		return *this;
	}
	for(int i = 0; i < this->getNTerms(); i++) {
		*this->terms[i] *= rh;
	}
	return *this;
}

GaussianSum &GaussianSum::operator/=(double rh) {
	for(int i = 0; i < this->getNTerms(); i++) {
		*this->terms[i] /= rh;
	}
	return *this;
}

GaussianSum &GaussianSum::operator-() const {
	GaussianSum *out = new GaussianSum(*this);
	for(int i = 0; i < this->getNTerms(); i++) {
		out->terms[i]->setCoef(-out->terms[i]->getCoef());
	}
	return *out;
}

GaussianSum &GaussianSum::operator+(const GaussianTerm &rh) const {
	GaussianSum *out = new GaussianSum(*this);
	*out += rh;

	return *out;
}

GaussianSum &GaussianSum::operator+(const GaussianSum &rh) const {
	GaussianSum *out = new GaussianSum(*this);
	*out += rh;

	return *out;
}

GaussianSum &GaussianSum::operator+(double rh) const {
	GaussianSum *out = new GaussianSum(*this);
	*out += rh;

	return *out;
}

GaussianSum &GaussianSum::operator-(const GaussianTerm &rh) const {
	GaussianSum *out = new GaussianSum(*this);
	*out -= rh;

	return *out;
}

GaussianSum &GaussianSum::operator-(const GaussianSum &rh) const {
	GaussianSum *out = new GaussianSum(*this);
	*out -= rh;

	return *out;
}

GaussianSum &GaussianSum::operator-(double rh) const {
	GaussianSum *out = new GaussianSum(*this);
	*out -= rh;

	return *out;
}

GaussianSum &GaussianSum::operator*(const GaussianTerm &rh) const {
	GaussianSum *out = new GaussianSum(*this);
	*out *= rh;

	return *out;
}

GaussianSum &GaussianSum::operator*(const GaussianSum &rh) const {
	GaussianSum *out = new GaussianSum(0);

	for(int i = 0; i < this->getNTerms(); i++) {
		for(int j = 0; j < rh.getNTerms(); j++) {
			*out += *this->terms[i] * *rh.terms[i];
		}
	}

	return *out;
}

GaussianSum &GaussianSum::operator*(double rh) const {
	GaussianSum *out = new GaussianSum(*this);
	*out *= rh;

	return *out;
}

GaussianSum &GaussianSum::operator/(double rh) const {
	GaussianSum *out = new GaussianSum(*this);
	*out /= rh;

	return *out;
}

GaussianTerm &operator*(double lh, const GaussianTerm &rh) {
	return rh * lh;
}

GaussianSum &operator+(double lh, const GaussianTerm &rh) {
	return GaussianSum(rh) + lh;
}

GaussianSum &operator+(double lh, const GaussianSum &rh) {
	return rh + lh;
}

GaussianSum &operator-(double lh, const GaussianTerm &rh) {
	return -rh + lh;
}

GaussianSum &operator-(double lh, const GaussianSum &rh) {
	return -rh + lh;
}

GaussianSum &operator*(double lh, const GaussianSum &rh) {
	return rh * lh;
}

static double choose(int n, int k) {
	return exp(lgamma(n + 1) - lgamma(n - k + 1) - lgamma(k + 1));
}

GaussianSum &translate(const GaussianTerm &term, const std::array<double, 3> &center) {
	GaussianSum *out = new GaussianSum();
	std::array<int, 3> expo;
	std::array<double, 3> newcent;
	double cx = term.getCenter(0) + center[0], cy = term.getCenter(1) + center[1], cz = term.getCenter(2) + center[2];

	newcent[0] = cx;
	newcent[1] = cy;
	newcent[2] = cz;

	for(int i = 0; i <= term.getExpo(0); i++) {
		expo[0] = i;
		double xcof = choose(term.getExpo(0), i) * pow(-center[0], term.getExpo(0) - i);
		if(xcof == 0) {
			continue;
		}
		for(int j = 0; j <= term.getExpo(1); j++) {
			expo[1] = j;
			double ycof = choose(term.getExpo(1), j) * pow(-center[1], term.getExpo(1) - j);
			if(ycof == 0) {
				continue;
			}
			for(int k = 0; k <= term.getExpo(2); k++) {
				double zcof = choose(term.getExpo(2), k) * pow(-center[2], term.getExpo(2) - k);
				if(zcof == 0) {
					continue;
				}

				expo[2] = k;

				*out += GaussianTerm(term.getCoef() * xcof * ycof * zcof, newcent, term.getMults(),
						expo);

			}
		}
	}

	return *out;
}

GaussianSum &compchem::translate(const GaussianSum &sum, const std::array<double, 3> &center) {
	GaussianSum *out = new GaussianSum();

	for(int i = 0; i < sum.getNTerms(); i++) {
		*out += translate(sum[i], center);
	}

	return *out;

}
