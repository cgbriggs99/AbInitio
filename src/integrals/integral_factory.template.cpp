#ifndef INTEGRAL_FACTORY_CPP
#define INTEGRAL_FACTORY_CPP


#include "integrals.hpp"
#include "../util/molecule.hpp"
#include "../util/tei_array.hpp"
#include "../util/atom.hpp"
#include <exception>
#include <stdexcept>
#include <vector>

template<typename Ints>
void compchem::IntegralFactory<Ints>::Smatrix(const compchem::Molecule *mol, double *out, int *dim) {
  if(out == nullptr || dim == nullptr) {
    throw new std::invalid_argument("The output values to the integral factory cn not be null.");
  }
  if(mol == nullptr) {
    throw new std::invalid_argument("The molecule input to the integral factory can not be null.");
  }
  int orbs = 0;
  for(int i = 0; i < mol->getsize(); i++) {
    orbs += mol->getatom(i).getnorbitals();
  }
  *dim = orbs;

  Ints method = Ints();
  int row = 0;
  for(int i = 0; i < mol->getsize(); i++) {
    for(int j = 0; j < mol->getatom(i).getnorbitals(); j++) {
      int col = 0;
      // fill in the diagonal.
      out[row * orbs + row] = method.overlap(
				static_cast<const compchem::GaussianOrbital *>(&mol->getatom(i).getorbital(j)),
				static_cast<const compchem::GaussianOrbital *>(&mol->getatom(i).getorbital(j)),
				mol->getatom(i).getpos(),
				mol->getatom(i).getpos());
      for(int k = 0; k <= i; k++) {
	if(k == i) {
	  for(int l = 0; l < j; l++) {
	    // Fill in the off-diagonal. Use the symmetry.
	    out[row * orbs + col] =
	      method.overlap(static_cast<const compchem::GaussianOrbital *>(&mol->getatom(i).getorbital(j)),
			     static_cast<const compchem::GaussianOrbital *>(&mol->getatom(k).getorbital(l)),
			     mol->getatom(i).getpos(),
			     mol->getatom(k).getpos());
	    out[col * orbs + row] = out[row * orbs + col];
	    col++;
	  }
	} else {
	  for(int l = 0; l < mol->getatom(k).getnorbitals(); l++) {
	    // Fill in the off-diagonal. Use the symmetry.
	    out[row * orbs + col] =
	      method.overlap(static_cast<const compchem::GaussianOrbital *>(&mol->getatom(i).getorbital(j)),
			     static_cast<const compchem::GaussianOrbital *>(&mol->getatom(k).getorbital(l)),
			     mol->getatom(i).getpos(),
			     mol->getatom(k).getpos());
	    out[col * orbs + row] = out[row * orbs + col];
	    col++;
	  }
	}
      }
      row++;
    }
  }
}
	  
			 
template<typename Ints>
void compchem::IntegralFactory<Ints>::Tmatrix(const compchem::Molecule *mol, double *out, int *dim) {
  if(out == nullptr || dim == nullptr) {
    throw new std::invalid_argument("The output values to the integral factory cn not be null.");
  }
  if(mol == nullptr) {
    throw new std::invalid_argument("The molecule input to the integral factory can not be null.");
  }
  int orbs = 0;
  for(int i = 0; i < mol->getsize(); i++) {
    orbs += mol->getatom(i).getnorbitals();
  }
  *dim = orbs;

  Ints method = Ints();
  int row = 0;
  for(int i = 0; i < mol->getsize(); i++) {
    for(int j = 0; j < mol->getatom(i).getnorbitals(); j++) {
      int col = 0;
      // fill in the diagonal.
      out[row * orbs + row] =
	method.laplacian(static_cast<const compchem::GaussianOrbital *>(&mol->getatom(i).getorbital(j)),
			 static_cast<const compchem::GaussianOrbital *>(&mol->getatom(i).getorbital(j)),
			 mol->getatom(i).getpos(),
			 mol->getatom(i).getpos());
      for(int k = 0; k <= i; k++) {
	if(k == i) {
	  for(int l = 0; l < j; l++) {
	    // Fill in the off-diagonal. Use the symmetry.
	    out[row * orbs + col] =
	      method.laplacian(static_cast<const compchem::GaussianOrbital *>(&mol->getatom(i).getorbital(j)),
			       static_cast<const compchem::GaussianOrbital *>(&mol->getatom(k).getorbital(l)),
			       mol->getatom(i).getpos(),
			       mol->getatom(k).getpos());
	    out[col * orbs + row] = out[row * orbs + col];
	    col++;
	  }
	} else {
	  for(int l = 0; l < mol->getatom(k).getnorbitals(); l++) {
	    // Fill in the off-diagonal. Use the symmetry.
	    out[row * orbs + col] =
	      method.laplacian(static_cast<const compchem::GaussianOrbital *>(&mol->getatom(i).getorbital(j)),
			       static_cast<const compchem::GaussianOrbital *>(&mol->getatom(k).getorbital(l)),
			       mol->getatom(i).getpos(),
			       mol->getatom(k).getpos());
	    out[col * orbs + row] = out[row * orbs + col];
	    col++;
	  }
	}
      }
      row++;
    }
  }
}

template<typename Ints>
void compchem::IntegralFactory<Ints>::Vmatrix(const compchem::Molecule *mol, double *out, int *dim) {
  if(out == nullptr || dim == nullptr) {
    throw new std::invalid_argument("The output values to the integral factory cn not be null.");
  }
  if(mol == nullptr) {
    throw new std::invalid_argument("The molecule input to the integral factory can not be null.");
  }
  int orbs = 0;
  for(int i = 0; i < mol->getsize(); i++) {
    orbs += mol->getatom(i).getnorbitals();
  }
  *dim = orbs;

  Ints method = Ints();
  int row = 0;
  for(int i = 0; i < mol->getsize(); i++) {
    for(int j = 0; j < mol->getatom(i).getnorbitals(); j++) {
      int col = 0;
      // fill in the diagonal.
      double sum = 0;
      for(int a = 0; a < mol->getsize(); a++) {
	sum += method.coulomb(static_cast<const compchem::GaussianOrbital *>(&mol->getatom(i).getorbital(j)),
			      static_cast<const compchem::GaussianOrbital *>(&mol->getatom(i).getorbital(j)),
			      mol->getatom(i).getpos(),
			      mol->getatom(i).getpos(),
			      mol->getatom(a));
      }
	
      out[row * orbs + row] = sum;
      for(int k = 0; k <= i; k++) {
	if(k == i) {
	  for(int l = 0; l < j; l++) {
	    // Fill in the off-diagonal. Use the symmetry.
	    sum = 0;
	    for(int a = 0; a < mol->getsize(); a++) {
	      sum += method.coulomb(static_cast<const compchem::GaussianOrbital *>(&mol->getatom(i).getorbital(j)),
				    static_cast<const compchem::GaussianOrbital *>(&mol->getatom(k).getorbital(l)),
				    mol->getatom(i).getpos(),
				    mol->getatom(k).getpos(),
				    mol->getatom(a));
	    }
	    out[col * orbs + row] = sum;
	    out[row * orbs + col] = sum;
	    col++;
	  }
	} else {
	  for(int l = 0; l < mol->getatom(k).getnorbitals(); l++) {
	    // Fill in the off-diagonal. Use the symmetry.
	    sum = 0;
	    for(int a = 0; a < mol->getsize(); a++) {
	      sum += method.coulomb(static_cast<const compchem::GaussianOrbital *>(&mol->getatom(i).getorbital(j)),
				    static_cast<const compchem::GaussianOrbital *>(&mol->getatom(k).getorbital(l)),
				    mol->getatom(i).getpos(),
				    mol->getatom(k).getpos(),
				    mol->getatom(a));
	    }
	    out[col * orbs + row] = sum;
	    out[row * orbs + col] = sum;
	    col++;
	  }
	}
      }
      row++;
    }
  }
}

template<typename Ints>
compchem::TEIArray *compchem::IntegralFactory<Ints>::TEIints(const compchem::Molecule *mol) {
  int orbs = 0;
  std::vector<const compchem::GaussianOrbital *> orbitals;
  std::vector<std::array<double, 3> > centers;
  
  for(int i = 0; i < mol->getsize(); i++) {
    orbs += mol->getatom(i).getnorbitals();
    for(int j = 0; j < mol->getatom(i).getnorbitals(); j++) {
      orbitals.push_back(static_cast<const compchem::GaussianOrbital *>
			 (&mol->getatom(i).getorbital(j)));
      centers.push_back(mol->getatom(i).getpos());
    }
  }
  compchem::TEIArray *out = new TEIArray(orbs);

  Ints method = Ints();

  int mu, nu, lam, sig;

  for(int i = 0; i < out->getsize(); i++) {
    out->indextoquad(i, &mu, &nu, &lam, &sig);

    out->at(mu, nu, lam, sig) = method.exchange(orbitals.at(mu),
						orbitals.at(nu),
						orbitals.at(lam),
						orbitals.at(sig),
						centers.at(mu),
						centers.at(nu),
						centers.at(lam),
						centers.at(sig));
  }

  orbitals.clear();
  centers.clear();
  
  return out;
}
  

#endif
