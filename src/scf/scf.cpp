#include "scf.hpp"

using namespace compchem;

SCFWfn::~SCFWfn() {
  if(S != nullptr) {
    delete[] S;
  }

  if(T != nullptr) {
    delete[] T;
  }

  if(V != nullptr) {
    delete[] V;
  }

  if(tei != nullptr) {
    delete tei;
  }
  
  if(Ca != nullptr) {
    delete[] Ca;
  }
  if(Cb != nullptr) {
    delete[] Cb;
  }
  if(Da != nullptr) {
    delete[] Da;
  }
  if(Db != nullptr) {
    delete[] Db;
  }
  if(Fa != nullptr) {
    delete[] Fa;
  }
  if(Fb != nullptr) {
    delete[] Fb;
  }
  if(es != nullptr) {
    delete[] es;
  }
}

const double *SCFWfn::getoverlap(int *dim) const {
  if(dim != nullptr) {
    *dim = this->dim;
  }
  return this->S;
}

const double *SCFWfn::getkinetic(int *dim) const {
  if(dim != nullptr) {
    *dim = this->dim;
  }
  return this->T;
}

const double *SCFWfn::getpotential(int *dim) const {
  if(dim != nullptr) {
    *dim = this->dim;
  }
  return this->V;
}

const TEIArray *SCFWfn::gettei() const {
  return this->tei;
}

void SCFWfn::setoverlap(double *data) {
  if(this->S != nullptr) {
    delete[] this->S;
  }

  this->S = data;
}

void SCFWfn::setkinetic(double *data) {
  if(this->T != nullptr) {
    delete[] this->T;
  }

  this->T = data;
}

void SCFWfn::setpotential(double *data) {
  if(this->V != nullptr) {
    delete[] this->V;
  }

  this->V = data;
}

void SCFWfn::settei(TEIArray *data) {
  if(this->tei != nullptr) {
    delete this->tei;
  }

  this->tei = data;
}


const double *SCFWfn::getcoefa(int *dim) const {
  if(dim != nullptr) {
    *dim = this->dim;
  }
  return this->Ca;
}

const double *SCFWfn::getcoefb(int *dim) const {
  if(dim != nullptr) {
    *dim = this->dim;
  }
  return this->Cb;
}

const double *SCFWfn::getdensa(int *dim) const {
  if(dim != nullptr) {
    *dim = this->dim;
  }
  return this->Da;
}

const double *SCFWfn::getdensb(int *dim) const {
  if(dim != nullptr) {
    *dim = this->dim;
  }
  return this->Db;
}

const double *SCFWfn::getfocka(int *dim) const {
  if(dim != nullptr) {
    *dim = this->dim;
  }
  return this->Fa;
}

const double *SCFWfn::getfockb(int *dim) const {
  if(dim != nullptr) {
    *dim = this->dim;
  }
  return this->Fb;
}

double SCFWfn::getenergy() const {
  return this->energy;
}

const double *SCFWfn::getenergies(int *dim) const {
  if(dim != nullptr) {
    *dim = this->dim;
  }
  return this->es;
}

int SCFWfn::getnorbs() const {
  return this->orbs;
}

void SCFWfn::setcoefa(double *arr) {
  if(this->Ca != nullptr) {
    delete[] this->Ca;
  }
  this->Ca = arr;
}

void SCFWfn::setcoefb(double *arr) {
  if(this->Cb != nullptr) {
    delete[] this->Cb;
  }
  this->Cb = arr;
}

void SCFWfn::setdensa(double *arr) {
  if(this->Da != nullptr) {
    delete[] this->Da;
  }
  this->Da = arr;
}

void SCFWfn::setdensb(double *arr) {
  if(this->Db != nullptr) {
    delete[] this->Db;
  }
  this->Db = arr;
}

void SCFWfn::setfocka(double *arr) {
  if(this->Fa != nullptr) {
    delete[] this->Fa;
  }
  this->Fa = arr;
}

void SCFWfn::setfockb(double *arr) {
  if(this->Fb != nullptr) {
    delete[] this->Fb;
  }
  this->Fb = arr;
}

void SCFWfn::setenergies(double *arr) {
  if(this->es != nullptr) {
    delete[] this->es;
  }
  this->es = arr;
}
