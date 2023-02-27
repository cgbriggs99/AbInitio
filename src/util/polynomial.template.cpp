
#ifndef __POLYNOMIAL_TCC__
#define __POLYNOMIAL_TCC__

#include <exception>
#include <vector>
#include <array>
#include "polynomial.hpp"
#include <cstdarg>
#include <cmath>
#include <cstring>
#include <stdexcept>
#include <cstdlib>
#include <initializer_list>

namespace compchem {

  // Done
template<int n>
Polynomial<n>::Polynomial(double scalar) {
  
  this->pows = (int *) std::calloc(n, sizeof(int));
  this->coefs = (double *) std::calloc(1, sizeof(double));
  this->size = 1;
  coefs[0] = scalar;
  for(int i = 0; i < n; i++) {
    this->pows[i] = 0;
  }
  this->reduce();
}

  // Done
template<int n>
Polynomial<n>::Polynomial(const std::vector<std::array<int, n> > &pows, const std::vector<double> &coefs) {
  this->pows = (int *) std::calloc(n * coefs.size(), sizeof(int));
  this->coefs = (double *) std::calloc(coefs.size(), sizeof(double));
  this->size = coefs.size();
  
  if(pows.size() != coefs.size()) {
    throw std::exception();
  }
  for(int i = 0; i < coefs.size(); i++) {
    this->coefs[i] = coefs[i];
    for(int j = 0; j < n; j++) {
      this->pows[i * n + j] = pows[i][j];
    }
  }
  this->reduce();
}

template<int n>
Polynomial<n>::Polynomial(const int *pows, const double *coefs, int size) {
  this->pows = (int *) std::calloc(n * size, sizeof(int));
  this->coefs = (double *) std::calloc(size, sizeof(double));
  this->size = size;

  std::memcpy((void *) this->pows, (void *) pows, n * size * sizeof(int));
  std::memcpy((void *) this->coefs, (void *) coefs, size * sizeof(double));

  this->reduce();
}

  // done
template<int n>
Polynomial<n>::Polynomial(const Polynomial<n> &copy) {
  this->pows = (int *) std::calloc(n * copy.getsize(), sizeof(int));
  this->coefs = (double *) std::calloc(copy.getsize(), sizeof(double));

  this->size = copy.getsize();

  if(copy.pows != nullptr) {
    std::memcpy((void *) this->pows, (void *) copy.pows,
		copy.getsize() * n * sizeof(int));
  }
  if(copy.coefs != nullptr) {
    std::memcpy((void *) this->coefs, (void *) copy.coefs,
		copy.getsize() * sizeof(double));
  }

  this->reduce();
}

  // done
template<int n>
Polynomial<n>::Polynomial() {
  this->pows = nullptr;
  this->coefs = nullptr;
  this->size = 0;
}

  // done
template<int n>
Polynomial<n>::~Polynomial() {
  if(this->coefs != nullptr) {
    std::free(this->coefs);
  }
  if(this->pows != nullptr) {
    std::free(this->pows);
  }
  this->coefs = nullptr;
  this->pows = nullptr;
}

  // done
template<int n>
const int *Polynomial<n>::gettermorder(int index) const {
  // This case is mathematically zero, but fails if let go.
  if(index == 0 && this->getsize() == 0) {
    return nullptr;
  }
  if(index < 0 || index >= this->getsize()) {
    throw new std::out_of_range("Array index out of range for polynomial power.");
  }
  return this->pows + index * n;
}

  // done
template<int n>
double Polynomial<n>::getcoef(int index) const {
  // This case is mathematically zero, but fails if let go.
  if(index == 0 && this->getsize() == 0) {
    return 0;
  }
  if(index < 0 || index >= this->getsize()) {
    throw new std::out_of_range("Array index out of range for polynomial coefficient.");
  }
  return this->coefs[index];
}

  // done
template<int n>
int Polynomial<n>::getsize() const {
  return this->size;
}

// done
template<int n>
void Polynomial<n>::reduce() {
  int end = this->size - 1;
  int temp_pow[n];
  // Filter the zero values to the end.
  for(int i = 0; i <= end; i++) {
    if(this->coefs[i] == 0 && i < end) {
      double temp = this->coefs[i];
      this->coefs[i] = this->coefs[end];
      this->coefs[end] = temp;
      
      std::memcpy((void *) temp_pow, (void *) (this->pows + i * n),
		  n * sizeof(int));
      std::memcpy((void *) (this->pows + i * n),
		  (void *) (this->pows + end * n), n * sizeof(int));
      std::memcpy((void *) (this->pows + end * n),
		  (void *) temp_pow, n * sizeof(int));
      i--;
      end--;
    } else if(this->coefs[i] == 0 && i == end) {
      end--;
    }
  }

  // Reallocate.
  if(end == -1) {
    std::free(this->pows);
    std::free(this->coefs);
    this->pows = nullptr;
    this->coefs = nullptr;
    this->size = 0;
  } else {
    void *ptr = std::realloc((void *) this->pows,
				      n * (end + 1) * sizeof(int));
    if(ptr == nullptr) {
      std::free(this->pows);
      std::free(this->coefs);
      this->pows = nullptr;
      this->coefs = nullptr;
      this->size = 0;
      return;
    } else {
      this->pows = (int *) ptr;
    }

    ptr = std::realloc((void *) this->coefs,
					  (end + 1) * sizeof(double));
    if(ptr == nullptr) {
      std::free(this->pows);
      std::free(this->coefs);
      this->pows = nullptr;
      this->coefs = nullptr;
      this->size = 0;
      return;
    } else {
      this->coefs = (double *) ptr;
    }
    this->size = end + 1;
  }
}

  // done
template<int n>
Polynomial<n> *Polynomial<n>::copy() const {
  return new Polynomial<n>(*this);
}

  // done
template<int n>
Polynomial<n> &Polynomial<n>::operator+=(double rh) {
  for(int i = 0; i < this->getsize(); i++) {
    int found = 1;
    for(int j = 0; j < n; j++) {
      if(this->pows[i * n + j] != 0) {
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
  this->pows = (int *) std::realloc((void *) this->pows,
			    n * (this->size + 1) * sizeof(int));
  this->coefs = (double *) std::realloc((void *) this->coefs,
			     (this->size + 1) * sizeof(double));
  for(int i = 0; i < n; i++) {
    this->pows[this->size * n + i] = 0;
  }
  this->coefs[this->size] = rh;
  this->size++;
  this->reduce();
  return *this;
}

  // done
template<int n>
Polynomial<n> &Polynomial<n>::operator-=(double rh) {
  for(int i = 0; i < this->getsize(); i++) {
    int found = 1;
    for(int j = 0; j < n; j++) {
      if(this->pows[i * n + j] != 0) {
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
  this->pows = (int *) std::realloc((void *) this->pows,
			    n * (this->size + 1) * sizeof(int));
  this->coefs = (double *) std::realloc((void *) this->coefs,
			     (this->size + 1) * sizeof(double));
  for(int i = 0; i < n; i++) {
    this->pows[this->size * n + i] = 0;
  }
  this->coefs[this->size] = -rh;
  this->size++;
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
  // Reallocate to limit the number of reallocations later.
  this->coefs = (double *) std::realloc((void *) this->coefs,
			     (this->getsize() + rh.getsize()) * sizeof(double));
  this->pows = (int *) std::realloc((void *) this->pows,
			    n * (this->getsize() + rh.getsize()) *
			    sizeof(double));
  int true_size = this->getsize();
  // loop through the right hand side
  for(int i = 0; i < rh.getsize(); i++) {
    int addto = 1;   // Whether the term has been added yet.
    for(int j = 0; j < this->getsize(); j++) {    // Loop through the left hand side
      int found = 1;    // Whether the left hand term is the same order as the right hand term.
      // Loop through the powers for the term.
      for(int k = 0; k < n; k++) {
	// If any power is not equal, the terms are different.
	if(this->pows[j * n + k] != rh.pows[i * n + k]) {
	  found = 0;
	  break;
	}
      }
      // Terms are the same order, so they can be added.
      if(found) {
	addto = 0;
	this->coefs[j] += rh.coefs[i];
	break;
      }
    }
    // Check if the term has been added. True if it has not been added.
    if(addto) {
      // Already reallocated the arrays. Copy data to the end.
      this->coefs[true_size] = rh.coefs[i];
      std::memcpy((void *) (this->pows + true_size * n),
		  (void *) (rh.pows + i * n), n * sizeof(int));
      true_size++;
    }
  }

  /*
   * By setting the size to the true_size without reallocating, the reduce
   * algorithm doesn't see past the true size. However, when the reduce
   * algorithm ends, it reallocates anyways to the stripped size, so there
   * should be no memory leaks (famous last words).
   */
  this->size = true_size;
  this->reduce();
  return *this;
}

template<int n>
Polynomial<n> &Polynomial<n>::operator-=(const Polynomial<n> &rh) {
  // Reallocate to limit the number of reallocations later.
  this->coefs = (double *) std::realloc((void *) this->coefs,
			     (this->getsize() + rh.getsize()) * sizeof(double));
  this->pows = (int *) std::realloc((void *) this->pows,
			    n * (this->getsize() + rh.getsize()) *
			    sizeof(double));
  int true_size = this->getsize();
  // loop through the right hand side
  for(int i = 0; i < rh.getsize(); i++) {
    int addto = 1;   // Whether the term has been added yet.
    for(int j = 0; j < this->getsize(); j++) {    // Loop through the left hand side
      int found = 1;    // Whether the left hand term is the same order as the right hand term.
      // Loop through the powers for the term.
      for(int k = 0; k < n; k++) {
	// If any power is not equal, the terms are different.
	if(this->pows[j * n + k] != rh.pows[i * n + k]) {
	  found = 0;
	  break;
	}
      }
      // Terms are the same order, so they can be subtracted.
      if(found) {
	addto = 0;
	this->coefs[j] -= rh.coefs[i];
	break;
      }
    }
    // Check if the term has been added. True if it has not been subtracted.
    if(addto) {
      // Already reallocated the arrays. Copy data to the end.
      this->coefs[true_size] = -rh.coefs[i];
      std::memcpy((void *) (this->pows + true_size * n),
		  (void *) (rh.pows + i * n), n * sizeof(int));
      true_size++;
    }
  }

  /*
   * By setting the size to the true_size without reallocating, the reduce
   * algorithm doesn't see past the true size. However, when the reduce
   * algorithm ends, it reallocates anyways to the stripped size, so there
   * should be no memory leaks (famous last words).
   */
  this->size = true_size;
  this->reduce();
  return *this;
}

template<int n>
Polynomial<n> &Polynomial<n>::operator*=(const Polynomial<n> &rh) {
  double *coefs = (double *) std::calloc(this->getsize() * rh.getsize(), sizeof(double));
  int *pows = (int *) std::calloc(n * this->getsize() * rh.getsize(), sizeof(int));
  int newsize = 0;

  // Cycle through left and right.
  for(int i = 0; i < this->getsize(); i++) {
    for(int j = 0; j < rh.getsize(); j++) {
      // Check to see if a term of the same power has been found.
      bool found = false;
      for(int k = 0; k < newsize; k++) {
	bool equal = true;
	for(int l = 0; l < n; l++) {
	  if(pows[k * n + l] != this->pows[i * n + l] + rh.pows[j * n + l]) {
	    equal = false;
	    break;
	  }
	}
	if(equal) {
	  found = true;
	  coefs[k] += this->coefs[i] * rh.coefs[j];
	  break;
	}
      }
      // Term has not been encountered. Create a new term.
      if(!found) {
	for(int k = 0; k < n; k++) {
	  pows[newsize * n + k] = this->pows[i * n + k] + rh.pows[j * n + k];
	}
	coefs[newsize] = this->coefs[i] * rh.coefs[j];
	newsize++;
      }
    }
  }

  std::free(this->coefs);
  std::free(this->pows);
  this->coefs = coefs;
  this->pows = pows;

  /*
   * By setting the size to the newsize without reallocating, the reduce
   * algorithm doesn't see past the true size. However, when the reduce
   * algorithm ends, it reallocates anyways to the stripped size, so there
   * should be no memory leaks (famous last words).
   */
  this->size = newsize;

  this->reduce();
  return *this;
}

template<int n>
Polynomial<n> Polynomial<n>::operator-() const {
  Polynomial<n> out = Polynomial(*this);
  for(int i = 0; i < out.getsize(); i++) {
    out.coefs[i] = -out.coefs[i];
  }
  out.reduce();
  return out;
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
      // constant term.
      if(this->pows[i * n + j] == 0) {
	prod *= 1;
      } else {
	prod *= std::pow(vars[j], this->pows[i * n + j]);
      }
    }
    sum += prod;
  }
  return sum;
}

template<int n>
double Polynomial<n>::eval(const double *r) const {
  double sum = 0;
  for(int i = 0; i < this->getsize(); i++) {
    double prod = this->coefs[i];
    for(int j = 0; j < n; j++) {
      if(this->pows[i * n + j] == 0) {
	prod *= 1;
      } else {
	prod *= std::pow(r[j], this->pows[i * n + j]);
      }
    }
    sum += prod;
  }
  return sum;
}

template<int n>
double Polynomial<n>::eval(std::initializer_list<double> r) const {
  double vars[n] = r;
  
  double sum = 0;
  for(int i = 0; i < this->getsize(); i++) {
    double prod = this->coefs[i];
    for(int j = 0; j < n; j++) {
      if(this->pows[i * n + j] == 0) {
	prod *= 1;
      } else {
	prod *= std::pow(vars[j], this->pows[i * n + j]);
      }
    }
    sum += prod;
  }
  return sum;
}

template<int m>
Polynomial<m> operator+(double lh, const Polynomial<m> &rh) {
  Polynomial<m> out = Polynomial<m>(rh);

  out += lh;
  
  return out;
}
  
template<int n>
Polynomial<n> operator+(Polynomial<n> lh, double rh) {
  Polynomial<n> out = lh;
  out += rh;
  return out;
}

template<int n>
Polynomial<n> operator+(Polynomial<n> lh, const Polynomial<n> &rh) {
  Polynomial<n> out = lh;
  out += rh;
  return out;
}

template<int m>
Polynomial<m> operator-(double lh, const Polynomial<m> &rh) {
  Polynomial<m> out = -rh;

  out += lh;
  
  return out;
}

template<int n>
Polynomial<n> operator-(Polynomial<n> lh, double rh) {
  Polynomial<n> out = lh;

  out -= rh;
  return out;
}

template<int n>
Polynomial<n> operator-(Polynomial<n> lh, const Polynomial<n> &rh) {
  Polynomial<n> out = lh;
  out -= rh;
  return out;
}

template<int m>
Polynomial<m> operator*(double lh, const Polynomial<m> &rh) {
  Polynomial<m> out = Polynomial<m>(rh);
  out *= lh;
  return out;
}

template<int n>
Polynomial<n> operator*(Polynomial<n> lh, double rh) {
  Polynomial<n> out = lh;
  out *= rh;
  return out;
}

template<int n>
Polynomial<n> operator*(Polynomial<n> lh, const Polynomial<n> &rh) {
  Polynomial<n> out = lh;
  out *= rh;
  return out;
}

template<int n>
Polynomial<n> operator/(Polynomial<n> lh, double rh) {
  Polynomial<n> out = lh;
  out /= rh;
  return out;
}

template<int n>
Polynomial<n> &Polynomial<n>::operator=(double d) {
  this->coefs = (double *) std::realloc((void *) this->coefs, sizeof(double));
  this->pows = (int *) std::realloc((void *) this->pows, n * sizeof(int));

  for(int i = 0; i < n; i++) {
    this->pows[i] = 0;
  }
  this->coefs[0] = d;
  this->size = 1;
  
  this->reduce();
  return *this;
}

template<int m>
Polynomial<m> pow(const Polynomial<m> &mant, int expo) {
  if(expo == 0 && mant != 0) {
    return Polynomial<m>(1);
  } else if(expo < 0) {
    throw new std::domain_error("Can not take a polynomial to a negative power at this time.");
  } else if(mant == 0 && expo == 0) {
    throw new std::domain_error("Can not raise the zero polynomial to the zero power.");
  } else if(mant == 0) {
    return Polynomial<m>(0);
  }
  
  Polynomial<m> out = Polynomial<m>(mant);
  // Start at 1, since we initialize to the mantissa.
  for(int i = 1; i < expo; i++) {
    out *= mant;
  }
  out.reduce();
  return out;
}

template<int n>
bool Polynomial<n>::operator==(double rh) const {
  if(this->getsize() > 1) {
    return false;
  } else if(this->getsize() == 0 && rh != 0) {
    return false;
  } else if(this->getsize() == 0 && rh == 0) {
    return true;
  } else if(this->getsize() == 1) {
    for(int i = 0; i < n; i++) {
      if(this->pows[i] != 0) {
        return false;
      }
    }
    return this->coefs[0] == rh;
  } else {
    return false;
  }
}

template<int n>
bool Polynomial<n>::operator!=(double rh) const {
  if(this->getsize() > 1) {
    return true;
  } else if(this->getsize() == 0 && rh != 0) {
    return true;
  } else if(this->getsize() == 0 && rh == 0) {
    return false;
  } else {
    for(int i = 0; i < n; i++) {
      if(this->pows[i] != 0) {
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
	  if(this->pows[i * m + k] != rh.pows[j * m + k]) {
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
  return !(*this == rh);
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
  double pos[n];

  va_list list;
  va_start(list, x);

  pos[0] = x;
  for(int i = 1; i < n; i++) {
    pos[i] = va_arg(list, double);
  }
  va_end(list);

  Polynomial<n> *out = new Polynomial<n>();
  Polynomial<n> *terms[n];
  // Set up a sneaky way to initialize the terms.
  int pow[2 * n - 1] = 0;
  pow[n - 1] = 1;
  double one = 1;
  for(int i = 0; i < n; i++) {
    terms[i] = new Polynomial<n>(pow + i, &one, 1);
    *terms[i] -= pos[i];
  }

  // Translate.
  for(int i = 0; i < this->getsize(); i++) {
    Polynomial<n> *term = new Polynomial<n>(this->coefs[i]);
    for(int j = 0; j < n; j++) {
      *term *= compchem::pow(*terms[i], this->pows[i * n + j]);
    }
    *out += *term;
    delete term;
  }
  // Clean up.
  delete this->coefs;
  this->coefs = out->coefs;
  out->coefs = nullptr;

  delete this->pows;
  this->pows = out->pows;
  out->pows = nullptr;

  this->size = out->size;
  delete out;
  delete[] terms;

  return *this;
    
}

template<int n>
Polynomial<n> &Polynomial<n>::operator=(const Polynomial<n> &rh) {
  if(this->coefs != nullptr) {
    std::free(this->coefs);
    this->coefs = nullptr;
  }
  if(this->pows != nullptr) {
    std::free(this->pows);
    this->pows = nullptr;
  }

  if(rh.coefs != nullptr) {
    this->coefs = (double *) std::calloc(rh.getsize(), sizeof(double));
    std::memcpy((void *) this->coefs, (void *) rh.coefs,
		rh.getsize() * sizeof(double));
  }
  if(rh.pows != nullptr) {
    this->pows = (int *) std::calloc(rh.getsize() * n, sizeof(int));
    std::memcpy((void *) this->pows, (void *) rh.pows,
		rh.getsize() * n *sizeof(int));
  }
  this->size = rh.getsize();
  return *this;
}

template<int n>
Polynomial<n> &Polynomial<n>::operator=(Polynomial<n> &&rh) {
  if(this->coefs != nullptr) {
    std::free(this->coefs);
    this->coefs = nullptr;
  }
  if(this->pows != nullptr) {
    std::free(this->pows);
    this->pows = nullptr;
  }

  if(rh.coefs != nullptr) {
    this->coefs = (double *) std::calloc(rh.getsize(), sizeof(double));
    std::memcpy((void *) this->coefs, (void *) rh.coefs,
		rh.getsize() * sizeof(double));
  }
  if(rh.pows != nullptr) {
    this->pows = (int *) std::calloc(rh.getsize() * n, sizeof(int));
    std::memcpy((void *) this->pows, (void *) rh.pows,
		rh.getsize() * n *sizeof(int));
  }
  this->size = rh.getsize();
  return *this;
}

template<int n>
bool Polynomial<n>::almost(double rh, double conv) const {
  if(this->getsize() > 1) {
    return false;
  } else if(this->getsize() == 0 && fabs(rh) > conv) {
    return false;
  } else if(this->getsize() == 0 && fabs(rh) <= conv) {
    return true;
  } else if(this->getsize() == 1) {
    for(int i = 0; i < n; i++) {
      if(this->pows[i] != 0) {
        return false;
      }
    }
    return fabs(this->coefs[0] - rh) <= conv;
  } else {
    return false;
  }
}

template<int n>
bool Polynomial<n>::almost(const Polynomial<n> &rh, double conv) const {
  if(this->getsize() != rh.getsize()) {
    return false;
  } else {
    for(int i = 0; i < this->getsize(); i++) {
      int found = 0;
      for(int j = 0; j < rh.getsize(); j++) {
        int equals = 1;
        for(int k = 0; k < n; k++) {
	  if(this->pows[i * n + k] != rh.pows[j * n + k]) {
	    equals = 0;
	    break;
	  }
	}
	if(equals == 0) {
	  continue;
	} else {
	  found = 1;
	  if(fabs(this->coefs[i] - rh.coefs[j]) > conv) {
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
void Polynomial<n>::setcoef(int index, double val) {
  if(index < 0 || index >= this->getsize()) {
    throw new std::out_of_range("Can not index polynomial.");
  }
  this->coefs[index] = val;
}
}
#endif
