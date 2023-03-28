#include "${CMAKE_SOURCE_DIR}/src/scf/rhf.hpp"
#include <cmath>
#include <${LAPACKE_HEADER}>
#include <${CBLAS_HEADER}>
#include <cstring>
#include <stdexcept>

using namespace compchem;
using namespace std;

const double *RHFWfn::getcoefa(int *dim) const {
  if(dim != nullptr) {
    *dim = this->dim;
  }
  return this->Ca;
}

const double *RHFWfn::getcoefb(int *dim) const {
  if(dim != nullptr) {
    *dim = this->dim;
  }
  return this->Ca;
}

const double *RHFWfn::getcoef(int *dim) const {
  if(dim != nullptr) {
    *dim = this->dim;
  }
  return this->Ca;
}

const double *RHFWfn::getdensa(int *dim) const {
  if(dim != nullptr) {
    *dim = this->dim;
  }
  return this->Da;
}

const double *RHFWfn::getdensb(int *dim) const {
  if(dim != nullptr) {
    *dim = this->dim;
  }
  return this->Da;
}

const double *RHFWfn::getdens(int *dim) const {
  if(dim != nullptr) {
    *dim = this->dim;
  }
  return this->Da;
}

const double *RHFWfn::getfocka(int *dim) const {
  if(dim != nullptr) {
    *dim = this->dim;
  }
  return this->Fa;
}

const double *RHFWfn::getfockb(int *dim) const {
  if(dim != nullptr) {
    *dim = this->dim;
  }
  return this->Fa;
}

const double *RHFWfn::getfock(int *dim) const {
  if(dim != nullptr) {
    *dim = this->dim;
  }
  return this->Fa;
}

void RHFWfn::setcoefa(double *arr) {
  if(this->Ca != nullptr) {
    delete[] this->Ca;
  }
  this->Ca = arr;
}

void RHFWfn::setcoefb(double *arr) {
  if(this->Ca != nullptr) {
    delete[] this->Ca;
  }
  this->Ca = arr;
}

void RHFWfn::setdensa(double *arr) {
  if(this->Da != nullptr) {
    delete[] this->Da;
  }
  this->Da = arr;
}

void RHFWfn::setdensb(double *arr) {
  if(this->Da != nullptr) {
    delete[] this->Da;
  }
  this->Da = arr;
}

void RHFWfn::setfocka(double *arr) {
  if(this->Fa != nullptr) {
    delete[] this->Fa;
  }
  this->Fa = arr;
}

void RHFWfn::setfockb(double *arr) {
  if(this->Fa != nullptr) {
    delete[] this->Fa;
  }
  this->Fa = arr;
}

void RHFWfn::setcoef(double *arr) {
  if(this->Ca != nullptr) {
    delete[] this->Ca;
  }
  this->Ca = arr;
}

void RHFWfn::setdens(double *arr) {
  if(this->Da != nullptr) {
    delete[] this->Da;
  }
  this->Da = arr;
}

void RHFWfn::setfock(double *arr) {
  if(this->Fa != nullptr) {
    delete[] this->Fa;
  }
  this->Fa = arr;
}

double RHF::energy(const Molecule *molecule,
		      const Wavefunction *wfn_in,
		      RHFWfn *out) const {

  const SCFWfn *wfn;

  try {
    wfn = static_cast<const SCFWfn *>(wfn_in);
  } catch(std::exception &e) {
    throw(e);
  }
  
  if(wfn->getelectrons() % 2 == 1 || wfn->getmultiplicity() != 1) {
    throw(*new std::invalid_argument("RHF can only handle singlet molecules."));
  }

  double enuc = nuclear_repulsion(*molecule);

  const double *S = wfn->getoverlap();
  double *eigreal = new double[wfn->getnorbs()];
  double *eigimag = new double[wfn->getnorbs()];
  double *eigbeta = new double[wfn->getnorbs()];
  double *energies = new double[wfn->getnorbs()];
  double *eigvecs = new double[wfn->getnorbs() * wfn->getnorbs()];
  int orbs = wfn->getnorbs(), elecs = wfn->getelectrons();
  int dim;
  int occ = wfn->getelectrons() / 2;
  double e0 = 0, e1 = 0;

  // Set up new matrices.
  double *D0 = new double[orbs * orbs], *D1 = new double[orbs * orbs];
  double *H = new double[orbs * orbs];
  double *C = new double[orbs * orbs];
  double *F = new double[orbs * orbs];

  // Create the hamiltonian and initial Fock matrix.
  const double *T = wfn->getkinetic();
  const double *V = wfn->getpotential();
  const TEIArray *tei = wfn->gettei();

  for(int i = 0; i < orbs * orbs; i++) {
    H[i] = T[i] + V[i];
    F[i] = T[i] + V[i];
  }

  double *work1 = new double[orbs * orbs],
    *work2 = new double[orbs * orbs];
  int cycles = 0,
    max_cycles = this->getoptions().getintoption("max scf cycles");

  double rms,
    rmsconv = this->getoptions().getdoubleoption("scf rms convergence"),
    enconv = this->getoptions().getdoubleoption("scf energy convergence");
  
  do {
    memcpy(work1, F, orbs * orbs * sizeof(double));
    memcpy(work2, S, orbs * orbs * sizeof(double));
    
    // Find the coefficients and energies.
    LAPACKE_dggev(LAPACK_ROW_MAJOR, 'N', 'V', orbs,
		  work1, orbs, work2, orbs, eigreal, eigimag, eigbeta,
		  nullptr, 1, eigvecs, orbs);
    
    // Sort the vectors.
    for(int i = 0; i < orbs; i++) {
      eigreal[i] /= eigbeta[i];
    }
    for(int i = 0; i < orbs; i++) {
      // Find the smallest not yet found.
      double min = eigreal[i];
      int minj = i;
      for(int j = i; j < orbs; j++) {
	if(eigreal[j] < min) {
	  minj = j;
	  min = eigreal[j];
	}
      }
      // Swap.
      energies[i] = eigreal[minj];
      double temp = eigreal[i];
      eigreal[i] = eigreal[minj];
      eigreal[minj] = temp;
      
      // Set the vector.
      for(int j = 0; j < orbs; j++) {
	C[j * orbs + i] = eigvecs[j * orbs + minj];
	temp = eigvecs[j * orbs + i];
	eigvecs[j * orbs + i] = eigvecs[j * orbs + minj];
	eigvecs[j * orbs + minj] = temp;
      }
    }

    // Compute the densities.
    for(int mu = 0; mu < orbs; mu++) {
      for(int nu = 0; nu <= mu; nu++) {
	D0[mu * orbs + nu] = D1[mu * orbs + nu];
	D1[mu * orbs + nu] = 0;
	for(int m = 0; m < occ; m++) {
	  D1[mu * orbs + nu] += C[mu * orbs + m] * C[nu * orbs + m];
	}
	D0[nu * orbs + mu] = D0[mu * orbs + nu];
	D1[nu * orbs + mu] = D1[mu * orbs + nu];
      }
    }

    // Compute the energy.
    e0 = e1;
    e1 = enuc;
    for(int mu = 0; mu < orbs; mu++) {
      for(int nu = 0; nu < orbs; nu++) {
	e1 += D1[mu * orbs + nu] * (H[mu * orbs + nu] + F[mu * orbs + nu]);
      }
    }

    // Compute the new Fock matrix.
    for(int mu = 0; mu < orbs; mu++) {
      for(int nu = 0; nu < orbs; nu++) {
	F[mu * orbs + nu] = H[mu * orbs + nu];
	for(int lam = 0; lam < orbs; lam++) {
	  for(int sig = 0; sig < orbs; sig++) {
	    F[mu * orbs + nu] += D1[lam * orbs + sig] *
	      (2 * tei->at(mu, nu, lam, sig) - tei->at(mu, lam, nu, sig));
	  }
	}
      }
    }

    // Compute the differences.
    rms = 0;
    for(int mu = 0; mu < orbs; mu++) {
      for(int nu = 0; nu < orbs; nu++) {
	rms += (D1[mu * orbs + nu] - D0[mu * orbs + nu]) *
	  (D1[mu * orbs + nu] - D0[mu * orbs + nu]);
      }
    }
    rms /= orbs * orbs;

    cycles++;
  } while(cycles < max_cycles &&
	  (cycles < 3 || rms > rmsconv || fabs(e1 - e0) > enconv));

  // Set the wavefunction return.
  if(out != nullptr) {
    out->setcoef(C);
    out->setdens(D1);
    out->setfock(F);
    out->setenergies(energies);
  } else {
    delete[] C;
    delete[] D1;
    delete[] F;
    delete[] energies;
  }
  delete[] eigreal;
  delete[] eigimag;
  delete[] eigbeta;
  delete[] eigvecs;
  delete[] D0;
  delete[] H;
  delete[] work1;
  delete[] work2;

  if(cycles >= max_cycles && (rms > rmsconv || fabs(e1 - e0) > enconv)) {
    throw(*new std::runtime_error("SCF energy did not converge!"));
  }

  return e1;
}
  
double RHF::energy(const Molecule *mol,
		      const Wavefunction *wfn_in) const {
  return this->energy(mol, wfn_in, nullptr);
}
