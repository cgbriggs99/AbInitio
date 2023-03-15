
#include "basis_set.hpp"
#include <cmath>
#include <vector>
#include <exception>
#include "../util/polynomial.hpp"
#include <cstdlib>
#include <cstring>
#include <stdexcept>

using namespace compchem;

GaussianOrbital::GaussianOrbital(int l, int ml, const std::vector<double> &coefs, const std::vector<double> &alphas) {
  this->l = l;
  this->ml = ml;
  this->coefs = new double[coefs.size()];
  this->alphas = new double[alphas.size()];
  this->size = coefs.size();
  this->norms = new double[coefs.size()];

  for(int i = 0; i < coefs.size(); i++) {
    this->coefs[i] = coefs[i];
    this->alphas[i] = alphas[i];
    this->norms[i] = 2 * std::pow(8 * alphas[i] * alphas[i] *
			      alphas[i] / M_PI, 0.25) *
      std::sqrt(2 * std::pow(8 * alphas[i], l) *
		std::exp(std::lgamma(l + 2) -
			 std::lgamma(2 * l + 3)));
  }
  this->harms = sphereharm(l, ml);
  if(coefs.size() != alphas.size()) {
    delete[] this->coefs;
    delete[] this->alphas;
    delete this->harms;
    delete[] this->norms;
    this->coefs = nullptr;
    this->alphas = nullptr;
    this->harms = nullptr;
    this->norms = nullptr;
    throw new std::exception();
  }
  this->sort();
}

GaussianOrbital::GaussianOrbital(const GaussianOrbital &copy) {
  this->l = copy.l;
  this->ml = copy.ml;
  this->coefs = new double[copy.getnterms()];
  this->alphas = new double[copy.getnterms()];
  this->norms = new double[copy.getnterms()];
  this->size = copy.getnterms();

  for(int i = 0; i < copy.getnterms(); i++) {
    this->coefs[i] = copy.coefs[i];
    this->alphas[i] = copy.alphas[i];
    this->norms[i] = copy.norms[i];
  }
  this->harms = copy.harms->copy();
}

GaussianOrbital::~GaussianOrbital() {
  if(this->coefs != nullptr) {
    delete[] this->coefs;
  }
  if(this->alphas != nullptr) {
    delete[] this->alphas;
  }
  if(this->harms != nullptr) {
    delete this->harms;
  }
  if(this->norms != nullptr) {
    delete this->norms;
  }
}

const double *GaussianOrbital::getcoefs() const {
  return this->coefs;
}

const double *GaussianOrbital::getalphas() const {
  return this->alphas;
}

double GaussianOrbital::getcoef(int index) const {
  if(index < 0 || index >= this->size) {
    throw new std::out_of_range("Term index out of bounds in Gaussian orbital.");
  }
  return this->coefs[index];
}

double GaussianOrbital::getalpha(int index) const {
  if(index < 0 || index >= this->size) {
    throw new std::out_of_range("Term index out of bounds in Gaussian orbital.");
  }
  return this->alphas[index];
}

int GaussianOrbital::getnterms() const {
  return this->size;
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
    sum += 2 * this->getcoef(i) * poly * exp(-this->getalpha(i) * r2) *
      std::pow(8 * this->getalpha(i) * this->getalpha(i) * this->getalpha(i)
	       / M_PI, 0.25) *
      std::sqrt(2 * std::pow(8 * this->getalpha(i), this->getl()) *
		std::exp(std::lgamma(this->getl() + 2) -
			 std::lgamma(2 * this->getl() + 3)));
  }
  return sum;
}

BasisOrbital *GaussianOrbital::copy() const {
  return new GaussianOrbital(*this);
}
  
GaussianOrbital &GaussianOrbital::operator=(const GaussianOrbital &other) {
  delete this->harms;
  this->harms = other.harms->copy();
  delete[] this->coefs;
  delete[] this->alphas;
  this->coefs = new double[other.getnterms()];
  this->alphas = new double[other.getnterms()];
  std::memcpy((void *) this->coefs, (void *) other.coefs,
	      other.getnterms() * sizeof(double));
  std::memcpy((void *) this->alphas, (void *) other.alphas,
	      other.getnterms() * sizeof(double));
  this->l = other.l;
  this->ml = other.ml;
  this->sort();
  return *this;
}

static void sort_internal(double *alphas, double *coefs, int size) {
  // Degenerate case.
  if(size <= 1) {
    return;
  } else if(size <= 16) {    // Small arrays get insertion sorted.
    for(int head = 0; head < size; head++) {
      double max = alphas[head];
      int maxi = head;
      for(int i = head; i < size; i++) {
	if(alphas[i] > max) {
	  max = alphas[i];
	  maxi = i;
	}
      }

      double temp = alphas[head];
      alphas[head] = alphas[maxi];
      alphas[maxi] = temp;
      temp = coefs[head];
      coefs[head] = coefs[maxi];
      coefs[maxi] = temp;
    }
  } else {    // Quick sort for large arrays.
    double pivot = 0;
    int split = 0;

    // Calculate the pivot.
    for(int i = 0; i < size; i++) {
      pivot += alphas[i];
    }
    pivot /= size;

    for(int i = 0; i < size; i++) {
      if(alphas[i] <= pivot) {
	if(i != split) {
	  double temp = alphas[split];
	  alphas[split] = alphas[i];
	  alphas[i] = temp;
	  temp = coefs[split];
	  coefs[split] = coefs[i];
	  coefs[i] = temp;
	}
	split++;
      }
    }
    // Sort the subarrays. Because it is done in-place, no need to copy back.
    sort_internal(alphas, coefs, split);
    sort_internal(alphas + split, coefs + split, size - split);
  }
}
    

void GaussianOrbital::sort(void) {
  // Sort by alphas.
  sort_internal(this->alphas, this->coefs, this->getnterms());
}
  
  
double GaussianOrbital::getnorm(int index) const {
  if(index < 0 || index >= this->size) {
    throw new std::out_of_range("Term index out of bounds in Gaussian orbital.");
  }
  return this->norms[index];
}

double GaussianOrbital::laplacian(double x, double y, double z) const {

  double sum = 0;

  for(int i = 0; i < this->getnterms(); i++) {
    for(int j = 0; j < this->getharms().getsize(); j++) {
      const int *pows1 = this->getharms().gettermorder(j);
      sum += this->getnorm(i) * this->getharms().getcoef(j) *
	this->getcoef(i);
	((pows1[0] >= 2? (pows1[0] * (pows1[0] - 1) *
			  std::pow(x, pows1[0] - 2) *
			  std::pow(y, pows1[1]) *
			  std::pow(z, pows1[2])) : 0) +
	 (pows1[1] >= 2? (pows1[1] * (pows1[1] - 1) *
			  std::pow(x, pows1[0]) *
			  std::pow(y, pows1[1] - 2) *
			  std::pow(z, pows1[2])) : 0) +
	 (pows1[2] >= 2? (pows1[2] * (pows1[2] - 1) *
			  std::pow(x, pows1[0]) *
			  std::pow(y, pows1[1]) *
			  std::pow(z, pows1[2] - 2)) : 0) -
	 4 * this->getalpha(i) * ((pows1[0] >= 1?
				 (pows1[0] * std::pow(x, pows1[0]) *
				  std::pow(y, pows1[1]) *
				  std::pow(z, pows1[2])): 0) +
				(pows1[1] >= 1?
				 (pows1[1] * std::pow(x, pows1[0]) *
				  std::pow(y, pows1[1]) *
				  std::pow(z, pows1[2])): 0) +
				(pows1[2] >= 1?
				 (pows1[2] * std::pow(x, pows1[0]) *
				  std::pow(y, pows1[1]) *
				  std::pow(z, pows1[2])): 0)) +
	 std::pow(x, pows1[0]) *
	 std::pow(y, pows1[1]) *
	 std::pow(z, pows1[2]) *
	 (4 * (x * x + y * y + z * z) * this->getalpha(i) * this->getalpha(i)
	  - 2 * this->getalpha(i))) *
	  std::exp(-this->getalpha(i) * (x * x + y * y + z * z));
    }
  }
  return sum;
}
