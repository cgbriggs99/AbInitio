
#include "basis_set.hpp"
#include <cmath>
#include <vector>
#include <exception>
#include "polynomial.hpp"

using namespace compchem;

GaussianOrbital::GaussianOrbital(int l, int ml, const std::vector<double> &coefs, const std::vector<double> &alphas) {
  this->l = l;
  this->ml = ml;
  this->coefs = std::vector<double>(coefs.size());
  this->alphas = std::vector<double>(alphas.size());

  for(int i = 0; i < coefs.size(); i++) {
    this->coefs[i] = coefs[i];
    this->alphas[i] = alphas[i];
  }
  this->harms = &sphereharm(l, ml);
  if(coefs.size() != alphas.size()) {
    delete this;
    throw new std::exception();
  }
}

GaussianOrbital::GaussianOrbital(const GaussianOrbital &copy) {
  this->l = copy.l;
  this->ml = copy.ml;
  this->coefs = std::vector<double>(copy.coefs.size());
  this->alphas = std::vector<double>(copy.alphas.size());

  for(int i = 0; i < coefs.size(); i++) {
    this->coefs[i] = copy.coefs[i];
    this->alphas[i] = copy.alphas[i];
  }
  this->harms = copy.harms;
}

const std::vector<double> &GaussianOrbital::getcoefs() const {
  return this->coefs;
}

const std::vector<double> &GaussianOrbital::getalphas() const {
  return this->alphas;
}

double GaussianOrbital::getcoef(int index) const {
  return this->coefs.at(index);
}

double GaussianOrbital::getalpha(int index) const {
  return this->alphas.at(index);
}

int GaussianOrbital::getnterms() const {
  return this->coefs.size();
}

int GaussianOrbital::getl() const {
  return this->l;
}

int GaussianOrbital::getml() const {
  return this->ml;
}

const Polynomial<3> &GaussianOrbital::getharms() const {
  return *this->harms;
}

double GaussianOrbital::eval(double x, double y, double z) const {
  double r2 = x * x + y * y + z * z;
  double poly = this->getharms().eval((double) x, (double) y, (double) z);

  double sum = 0;
  for(int i = 0; i < this->getnterms(); i++) {
    sum += this->getcoef(i) * std::pow(2 * this->getalpha(i) / M_PI, 0.75) * poly * exp(-this->getalpha(i) * r2);
  }
  return sum;
}

BasisOrbital &GaussianOrbital::copy() const {
  return *new GaussianOrbital(*this);
}
  
