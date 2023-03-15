
#include "basis_set.hpp"
#include <cmath>
#include "../util/polynomial.hpp"



using namespace compchem;

double compchem::slater_rule(int n, int l, const compchem::GSConfig &conf) {
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

SlaterOrbital::SlaterOrbital(double Zeff, int n, int l, int ml) :
  Zeff(Zeff), n(n), l(l), ml(ml), harms(sphereharm(l, ml)) {
  ;
}

SlaterOrbital::SlaterOrbital(const SlaterOrbital &copy) :
  Zeff(copy.Zeff), n(copy.n), l(copy.l), ml(copy.ml), harms(copy.harms->copy()) {
  ;
}

SlaterOrbital::~SlaterOrbital() {
  delete this->harms;
}

double SlaterOrbital::eval(double x, double y, double z) const {

  double r = std::sqrt(x * x + y * y + z * z);

  double harm = this->harms->eval(x, y, z);
  if(this->getn() - this->getl() <= 1) {
    return std::exp(-r * this->getZeff()) * harm *
      std::pow(2 * this->getZeff(), this->getn()) *
      std::sqrt(2 * this->getZeff() / std::tgamma(2 * this->getn() + 1));
  }

  return std::pow(r, this->getn() - this->getl() - 1) *
    std::exp(-this->getZeff() * r) * harm *
    std::pow(2 * this->getZeff(), this->getn()) *
    std::sqrt(2 * this->getZeff() / std::tgamma(2 * this->getn() + 1));
}


  
double SlaterOrbital::getZeff() const {
  return this->Zeff;
}

int SlaterOrbital::getn() const {
  return this->n;
}

int SlaterOrbital::getl() const {
  return this->l;
}

int SlaterOrbital::getml() const {
  return this->getml();
}

const Polynomial<3> &SlaterOrbital::getharms() const {
  return *this->harms;
}

BasisOrbital *SlaterOrbital::copy() const {
  return new SlaterOrbital(*this);
}

double SlaterOrbital::laplacian(double x, double y, double z) const {
  double sum = 0;
  for(int i = 0; i < this->getharms().getsize(); i++) {
    const int *pows1 = this->getharms().gettermorder(i);
    sum += ((pows1[0] >= 2? (pows1[0] * (pows1[0] - 1) *
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
	    2 * this->getZeff() / std::sqrt(x * x + y * y + z * z) *
	    ((pows1[0] >= 1? (pows1[0] * std::pow(x, pows1[0])): 0) +
	     (pows1[1] >= 1? (pows1[1] * std::pow(y, pows1[1])): 0) +
	     (pows1[2] >= 1? (pows1[2] * std::pow(z, pows1[2])): 0)) +
	    0.5 * std::pow(x, pows1[0]) *
	    std::pow(y, pows1[1]) *
	    std::pow(z, pows1[2]) * (this->getZeff() * this->getZeff() *
				     (x * x + y * y + z * z) -
				     2.5 * this->getZeff()) /
	    std::sqrt(x * x + y * y + z * z)) *
      std::exp(-this->getZeff() * std::sqrt(x * x + y * y + z * z));
  }
  return std::pow(2 * this->getZeff(), this->getn()) *
    std::sqrt(2 * this->getZeff() / std::tgamma(2 * this->getn() + 1)) *
    sum;
}
