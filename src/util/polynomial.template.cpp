
#ifndef __POLYNOMIAL_TCC__
#define __POLYNOMIAL_TCC__

#include <exception>
#include <vector>
#include <array>
#include "polynomial.hpp"
#include <cstdarg>
#include <cmath>

namespace compchem {

template<int n>
Polynomial<n>::Polynomial(double scalar) {
  this->pows = std::vector<std::array<int, n> >(1);
  this->coefs = std::vector<double>(1);

  std::array<int, n> zero;
  for(int i = 0; i < n; i++) {
    zero[i] = 0;
  }

  this->pows.push_back(zero);
  this->coefs.push_back(scalar);
}

template<int n>
Polynomial<n>::Polynomial(const std::vector<std::array<int, n> > &pows, const std::vector<double> &coefs) {
  if(pows.size() != coefs.size()) {
    throw std::exception();
  }
  this->pows = std::vector<std::array<int, n> >(pows.size());
  this->coefs = std::vector<double>(coefs.size());

  for(int i = 0; i < coefs.size(); i++) {
    this->pows[i] = pows[i];
    this->coefs[i] = coefs[i];
  }
}

template<int n>
Polynomial<n>::Polynomial(const Polynomial<n> &copy) {
  this->pows = std::vector<std::array<int, n> >(copy.pows.size());
  this->coefs = std::vector<double>(copy.coefs.size());

  for(int i = 0; i < copy.coefs.size(); i++) {
    this->pows[i] = copy.pows[i];
    this->coefs[i] = copy.coefs[i];
  }
}

template<int n>
Polynomial<n>::Polynomial() {
  this->pows = std::vector<std::array<int, n> >(0);
  this->coefs = std::vector<double>(0);
}

template<int n>
const std::array<int, n> &Polynomial<n>::gettermorder(int index) const {
  return this->pows.at(index);
}

template<int n>
double Polynomial<n>::getcoef(int index) const {
  return this->coefs.at(index);
}

template<int n>
int Polynomial<n>::getsize() const {
  return this->coefs.size();
}

template<int n>
void Polynomial<n>::reduce() {
  for(int i = 0; i < this->coefs.size(); i++) {
    if(this->coefs[i] == 0) {
      this->coefs.erase(this->coefs.begin() + i);
      this->pows.erase(this->pows.begin() + i);
      i--;
    }
  }
}

template<int n>
Polynomial<n> *Polynomial<n>::copy() const {
  return new Polynomial<n>(*this);
}

template<int n>
Polynomial<n> &Polynomial<n>::operator+=(double rh) {
  for(int i = 0; i < this->getsize(); i++) {
    int found = 1;
    for(int j = 0; j < n; j++) {
      if(this->pows[i][j] != 0) {
	found = 0;
	break;
      }
    }
    if(found) {
      this->coefs[i] += rh;
      this->reduce();
      return *this;
    }
  }
  std::array<int, n> zero;
  for(int i = 0; i < n; i++) {
    zero[i] = 0;
  }
  this->pows.push_back(zero);
  this->coefs.push_back(rh);
  this->reduce();
  return *this;
}

template<int n>
Polynomial<n> &Polynomial<n>::operator-=(double rh) {
  for(int i = 0; i < this->getsize(); i++) {
    int found = 1;
    for(int j = 0; j < n; j++) {
      if(this->pows[i][j] != 0) {
	found = 0;
	break;
      }
    }
    if(found) {
      this->coefs[i] -= rh;
      this->reduce();
      return *this;
    }
  }
  std::array<int, n> zero;
  for(int i = 0; i < n; i++) {
    zero[i] = 0;
  }
  this->pows.push_back(zero);
  this->coefs.push_back(-rh);
  this->reduce();
  return *this;
}

template<int n>
Polynomial<n> &Polynomial<n>::operator*=(double rh) {
  for(int i = 0; i < this->getsize(); i++) {
    this->coefs[i] *= rh;
  }
  this->reduce();
  return *this;
}

template<int n>
Polynomial<n> &Polynomial<n>::operator/=(double rh) {
  for(int i = 0; i < this->getsize(); i++) {
    this->coefs[i] /= rh;
  }
  this->reduce();
  return *this;
}

template<int n>
Polynomial<n> &Polynomial<n>::operator+=(const Polynomial<n> &rh) {
  for(int i = 0; i < rh.getsize(); i++) {
    int addto = 1;
    for(int j = 0; j < this->getsize(); j++) {
      int found = 1;
      for(int k = 0; k < n; k++) {
	if(this->pows[j][k] != rh.pows[i][k]) {
	  found = 0;
	  break;
	}
      }
      if(found) {
	addto = 0;
	this->coefs[j] += rh.coefs[i];
	break;
      }
    }
    if(addto) {
      this->coefs.push_back(rh.coefs[i]);
      this->pows.push_back(rh.pows[i]);
    }
  }
  this->reduce();
  return *this;
}

template<int n>
Polynomial<n> &Polynomial<n>::operator-=(const Polynomial<n> &rh) {
  for(int i = 0; i < rh.getsize(); i++) {
    int addto = 1;
    for(int j = 0; j < this->getsize(); j++) {
      int found = 1;
      for(int k = 0; k < n; k++) {
	if(this->pows[j][k] != rh.pows[i][k]) {
	  found = 0;
	  break;
	}
      }
      if(found) {
	addto = 0;
	this->coefs[j] -= rh.coefs[i];
	break;
      }
    }
    if(addto) {
      this->coefs.push_back(-rh.coefs[i]);
      this->pows.push_back(rh.pows[i]);
    }
  }
  this->reduce();
  return *this;
}

template<int n>
Polynomial<n> &Polynomial<n>::operator*=(const Polynomial<n> &rh) {
  Polynomial<n> *poly = new Polynomial<n>();
  for(int i = 0; i < this->getsize(); i++) {
    for(int j = 0; j < rh.getsize(); j++) {
      std::array<int, n> power;

      for(int k = 0; k < n; k++) {
	power[k] = this->pows[i][k] + rh.pows[j][k];
      }
      poly->pows.push_back(power);
      poly->coefs.push_back(this->coefs[i] * rh.coefs[j]);
    }
  }
  this->coefs.swap(poly->coefs);
  this->pows.swap(poly->pows);
  delete poly;

  this->reduce();
  return *this;
}

template<int n>
Polynomial<n> &Polynomial<n>::operator-() const {
  Polynomial<n> *out = new Polynomial(*this);
  for(int i = 0; i < out->getsize(); i++) {
    out->coefs[i] = -out->coefs[i];
  }
  out->reduce();
  return *out;
}

template<int n>
double Polynomial<n>::eval(double x,...) const {
  va_list list;
  va_start(list, x);

  std::array<double, n> vars;
  vars[0] = x;
  for(int i = 1; i < n; i++) {
    vars[i] = va_arg(list, double);
  }
  va_end(list);

  double sum = 0;

  for(int i = 0; i < this->getsize(); i++) {
    double prod = this->coefs[i];
    for(int j = 0; j < n; j++) {
      prod *= std::pow(vars[j], this->pows[i][j]);
    }
    sum += prod;
  }
  return sum;
}

template<int m>
Polynomial<m> &operator+(double lh, const Polynomial<m> &rh) {
  Polynomial<m> *out = new Polynomial<m>(rh);

  *out += lh;
  return *out;
}
  
template<int m>
Polynomial<m> &operator+(Polynomial<m> lh, double rh) {
  Polynomial<m> *out = new Polynomial<m>(lh);

  *out += rh;
  return *out;
}

template<int m>
Polynomial<m> &operator+(Polynomial<m> lh, const Polynomial<m> &rh) {
  Polynomial<m> *out = new Polynomial<m>(rh);
  *out += lh;
  return *out;
}

template<int m>
Polynomial<m> &operator-(double lh, const Polynomial<m> &rh) {
  Polynomial<m> *out = new Polynomial<m>(-rh);
  *out += lh;
  return *out;
}

template<int m>
Polynomial<m> &operator-(Polynomial<m> lh, double rh) {
  Polynomial<m> *out = new Polynomial<m>(lh);

  *out -= rh;
  return out;
}

template<int m>
Polynomial<m> &operator-(Polynomial<m> lh, const Polynomial<m> &rh) {
  Polynomial<m> *out = new Polynomial<m>(lh);
  *out -= rh;
  return *out;
}

template<int m>
Polynomial<m> &operator*(double lh, const Polynomial<m> &rh) {
  Polynomial<m> *out = new Polynomial<m>(rh);
  *out *= lh;
  return *out;
}

template<int m>
Polynomial<m> &operator*(Polynomial<m> lh, double rh) {
  Polynomial<m> *out = new Polynomial<m>(lh);
  *out *= rh;
  return *out;
}

template<int m>
Polynomial<m> &operator*(Polynomial<m> lh, const Polynomial<m> &rh) {
  Polynomial<m> *out = new Polynomial<m>(rh);
  *out *= lh;
  return *out;
}

template<int m>
Polynomial<m> &operator/(Polynomial<m> lh, double rh) {
  Polynomial<m> *out = new Polynomial<m>(lh);
  *out *= rh;
  return *out;
}

template<int n>
Polynomial<n> &Polynomial<n>::operator=(double d) {
  this->coefs.empty();
  this->pows.empty();

  std::array<int, n> zero;
  for(int i = 0; i < n; i++) {
    zero[i] = 0;
  }
  this->coefs.push_back(d);
  this->pows.push_back(zero);
  this->reduce();
  return *this;
}

template<int m>
Polynomial<m> &pow(const Polynomial<m> &mant, int expo) {
  if(expo == 0 && mant != 0) {
    return *new Polynomial<m>(1);
  } else if(expo < 0 || mant == 0) {
    throw new std::exception();
  }
  
  Polynomial<m> *out = new Polynomial<m>(mant);
  for(int i = 0; i < expo; i++) {
    *out *= mant;
  }
  out->reduce();
  return *out;
}

template<int n>
bool Polynomial<n>::operator==(double rh) const {
  if(this->getsize() > 1) {
    return false;
  } else if(this->getsize() == 0 && rh != 0) {
    return false;
  } else {
    for(int i = 0; i < n; i++) {
      if(this->pows[0][i] != 0) {
        return false;
      }
    }
    return this->coefs[0] == rh;
  }
}

template<int n>
bool Polynomial<n>::operator!=(double rh) const {
  if(this->getsize() > 1) {
    return true;
  } else if(this->getsize() == 0 && rh != 0) {
    return true;
  } else {
    for(int i = 0; i < n; i++) {
      if(this->pows[0][i] != 0) {
        return true;
      }
    }
    return this->coefs[0] != rh;
  }
}

template<int m>
bool Polynomial<m>::operator==(const Polynomial<m> &rh) const {
  if(this->getsize() != rh.getsize()) {
    return false;
  } else {
    for(int i = 0; i < this->getsize(); i++) {
      int found = 0;
      for(int j = 0; j < rh.getsize(); j++) {
        int equals = 1;
        for(int k = 0; k < m; k++) {
	  if(this->pows[i][k] != rh.pows[j][k]) {
	    equals = 0;
	    break;
	  }
	}
	if(equals == 0) {
	  continue;
	} else {
	  found = 1;
	  if(this->coefs[i] != rh.coefs[j]) {
	    return false;
	  } else {
	    break;
	  }
	}
	if(found == 0) {
	  return false;
	}
      }
    }
    return true;
  }
}

template<int n>
bool Polynomial<n>::operator!=(const Polynomial &rh) const {
  return *this != rh;
}
	

template<int m>
bool operator==(double lh, const Polynomial<m> &rh) {
  return rh == lh;
}

template<int m>
bool operator!=(double lh, const Polynomial<m> &rh) {
  return rh != lh;
}

template<int n>
Polynomial<n> &Polynomial<n>::translate(double x, ...) {
  std::array<double, n> pos;

  va_list list;
  va_start(list, x);

  pos[0] = x;
  for(int i = 1; i < n; i++) {
    pos[i] = va_arg(list, double);
  }
  va_end(list);

  std::vector<std::array<int, n> > pow = std::vector<std::array<int, n> >(2);
  std::vector<double> coefs = std::vector<double>(2);

  std::vector<Polynomial<n> *> terms = std::vector<Polynomial<n> *>(n);

  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      if(i == j) {
        pow[1][j] = 1;
      } else {
        pow[1][j] = 0;
      }
      pow[0][j] = 0;
    }
    coefs[1] = 1;
    coefs[0] = -pos[i];
    terms.push_back(new Polynomial<n>(pow, coefs));
  }

  Polynomial<n> *out = new Polynomial<n>(),
    *polyterm = new Polynomial<n>(1);

  for(int i = 0; i < this->getsize(); i++) {
    *polyterm = this->getcoef(i);

    for(int j = 0; j < n; j++) {
      *polyterm *= pow(terms[j], this->getpow(i)[j]);
    }
    *out += *polyterm;
  }

  delete polyterm;
  for(int i = 0; i < n; i++) {
    delete terms[i];
  }
  this->coefs.swap(out->coefs);
  this->pows.swap(out->pows);
  delete out;
  return *this;
}

}
#endif
