#include "tei_array.hpp"
#include <exception>
#include <stdexcept>
#include <cmath>

using namespace compchem;

TEIArray::TEIArray(double * __restrict__ data, int dim) {
  this->data = data;
  this->dim = dim;
}

TEIArray::TEIArray(int dim) {
  int mn = (dim * (dim - 1)) / 2 + dim - 1,
    mnls = (mn * (mn + 1)) / 2 + mn;
  this->data = new double[mnls + 1];
  this->dim = dim;
}

TEIArray::TEIArray(const TEIArray &copy) {
  this->dim = copy.getdim();
  int mn = (this->dim * (this->dim + 1)) / 2 + this->dim,
    mnls = (mn * (mn + 1)) / 2 + mnls;
  this->data = new double[mnls];
}

TEIArray::~TEIArray() {
  delete[] this->data;
}

double TEIArray::at(int mu, int nu, int lam, int sig) const {
  if(mu < 0 || nu < 0 || lam < 0 || sig < 0 ||
     mu >= this->dim || nu >= this->dim ||
     lam >= this->dim || sig >= this->dim) {
    throw *new std::out_of_range("Array index out of range for two-electron integrals.");
  }
  int m, n, l, s, mn, ls, mnls;
  if(mu >= nu) {
    m = mu;
    n = nu;
  } else {
    m = nu;
    n = mu;
  }
  if(lam >= sig) {
    l = lam;
    s = sig;
  } else {
    l = sig;
    s = lam;
  }
  mn = (m * (m + 1)) / 2 + n;
  ls = (l * (l + 1)) / 2 + s;

  if(mn < ls) {
    int temp = mn;
    mn = ls;
    ls = temp;
  }
  mnls = (mn * (mn + 1)) / 2 + ls;
  return this->data[mnls];
}

double &TEIArray::at(int mu, int nu, int lam, int sig) {
  if(mu < 0 || nu < 0 || lam < 0 || sig < 0 ||
     mu >= this->dim || nu >= this->dim ||
     lam >= this->dim || sig >= this->dim) {
    throw *new std::out_of_range("Array index out of range for two-electron integrals.");
  }
  int m, n, l, s, mn, ls, mnls;
  if(mu >= nu) {
    m = mu;
    n = nu;
  } else {
    m = nu;
    n = mu;
  }
  if(lam >= sig) {
    l = lam;
    s = sig;
  } else {
    l = sig;
    s = lam;
  }
  mn = (m * (m + 1)) / 2 + n;
  ls = (l * (l + 1)) / 2 + s;

  if(mn < ls) {
    int temp = mn;
    mn = ls;
    ls = temp;
  }
  mnls = (mn * (mn + 1)) / 2 + ls;
  return this->data[mnls];
}

double TEIArray::at_direct(int index) const {
  if(index < 0 || index > this->getsize()) {
    throw *new std::out_of_range("Array index out of range for two-electron integrals.");
  }

  return this->data[index];
}

double &TEIArray::at_direct(int index) {
  if(index < 0 || index > this->getsize()) {
    throw *new std::out_of_range("Array index out of range for two-electron integrals.");
  }

  return this->data[index];
}

double TEIArray::operator()(int mu, int nu, int lam, int sig) const {
  if(mu < 0 || nu < 0 || lam < 0 || sig < 0 ||
     mu >= this->dim || nu >= this->dim ||
     lam >= this->dim || sig >= this->dim) {
    throw new std::out_of_range("Array index out of range for two-electron integrals.");
  }
  int m, n, l, s, mn, ls, mnls;
  if(mu >= nu) {
    m = mu;
    n = nu;
  } else {
    m = nu;
    n = mu;
  }
  if(lam >= sig) {
    l = lam;
    s = sig;
  } else {
    l = sig;
    s = lam;
  }
  mn = (m * (m + 1)) / 2 + n;
  ls = (l * (l + 1)) / 2 + s;

  if(mn < ls) {
    int temp = mn;
    mn = ls;
    ls = temp;
  }
  mnls = (mn * (mn + 1)) / 2 + ls;
  return this->data[mnls];
}

double &TEIArray::operator()(int mu, int nu, int lam, int sig) {
  if(mu < 0 || nu < 0 || lam < 0 || sig < 0 ||
     mu >= this->dim || nu >= this->dim ||
     lam >= this->dim || sig >= this->dim) {
    throw new std::out_of_range("Array index out of range for two-electron integrals.");
  }
  int m, n, l, s, mn, ls, mnls;
  if(mu >= nu) {
    m = mu;
    n = nu;
  } else {
    m = nu;
    n = mu;
  }
  if(lam >= sig) {
    l = lam;
    s = sig;
  } else {
    l = sig;
    s = lam;
  }
  mn = (m * (m + 1)) / 2 + n;
  ls = (l * (l + 1)) / 2 + s;

  if(mn < ls) {
    int temp = mn;
    mn = ls;
    ls = temp;
  }
  mnls = (mn * (mn + 1)) / 2 + ls;
  return this->data[mnls];
}

const double *TEIArray::getdata() const {
  return this->data;
}

int TEIArray::getdim() const {
  return this->dim;
}

int TEIArray::getindex(int mu, int nu, int lam, int sig) const {
  if(mu < 0 || nu < 0 || lam < 0 || sig < 0 ||
     mu >= this->dim || nu >= this->dim ||
     lam >= this->dim || sig >= this->dim) {
    throw new std::out_of_range("Array index out of range for two-electron integrals.");
  }
  int m, n, l, s, mn, ls;
  if(mu >= nu) {
    m = mu;
    n = nu;
  } else {
    m = nu;
    n = mu;
  }
  if(lam >= sig) {
    l = lam;
    s = sig;
  } else {
    l = sig;
    s = lam;
  }
  mn = (m * (m + 1)) / 2 + n;
  ls = (l * (l + 1)) / 2 + s;

  if(mn < ls) {
    int temp = mn;
    mn = ls;
    ls = temp;
  }
  return (mn * (mn + 1)) / 2 + ls;
}

int TEIArray::getsize() const {
  int mn = (this->getdim() * (this->getdim() - 1)) / 2 + this->getdim() - 1;

  return (mn * (mn + 1)) / 2 + mn + 1;
}

void TEIArray::indextoquad(int index, int *mu, int *nu, int *lam, int *sig) const {
  if(index < 0 || index >= this->getsize()) {
    throw new std::out_of_range("Array index out of range for two-electron integrals.");
  }
  // Find the largest triangular number that is less than or equal to
  // the index. This corresponds to the mu,nu index.
  int mn = (int) std::floor((std::sqrt(1.0 + 8.0 * index) - 1.0) / 2.0);
  int Tmn = (mn * (mn + 1)) / 2;

  // The lam,sig index is just the difference between the index and this
  // triangular number.
  int ls = index - Tmn;

  // Find the mu and nu indices in a similar way.
  *mu = (int) std::floor((std::sqrt(1.0 + 8.0 * mn) - 1.0) / 2.0);
  int Tmu = (*mu * (*mu + 1)) / 2;
  *nu = mn - Tmu;

  // Find the lam and sig indices in a similar way.
  *lam = (int) std::floor((std::sqrt(1.0 + 8.0 * ls) - 1.0) / 2.0);
  int Tlam = (*lam * (*lam + 1)) / 2;
  *sig = ls - Tlam;
}
  
