#include "${CMAKE_SOURCE_DIR}/src/scf/uhf.hpp"
#include <cmath>
#include <${LAPACKE_HEADER}>
#include <${CBLAS_HEADER}>
#include <cstring>
#include <stdexcept>
#include <cstdio>

using namespace compchem;
using namespace std;

const double *UHFWfn::getcoefa(int *dim) const {
  if(dim != nullptr) {
    *dim = this->dim;
  }
  return this->Ca;
}

const double *UHFWfn::getcoefb(int *dim) const {
  if(dim != nullptr) {
    *dim = this->dim;
  }
  return this->Cb;
}

const double *UHFWfn::getdensa(int *dim) const {
  if(dim != nullptr) {
    *dim = this->dim;
  }
  return this->Da;
}

const double *UHFWfn::getdensb(int *dim) const {
  if(dim != nullptr) {
    *dim = this->dim;
  }
  return this->Db;
}

const double *UHFWfn::getfocka(int *dim) const {
  if(dim != nullptr) {
    *dim = this->dim;
  }
  return this->Fa;
}

const double *UHFWfn::getfockb(int *dim) const {
  if(dim != nullptr) {
    *dim = this->dim;
  }
  return this->Fb;
}

void UHFWfn::setcoefa(double *arr) {
  if(this->Ca != nullptr) {
    delete[] this->Ca;
  }
  this->Ca = arr;
}

void UHFWfn::setcoefb(double *arr) {
  if(this->Cb != nullptr) {
    delete[] this->Cb;
  }
  this->Cb = arr;
}

void UHFWfn::setdensa(double *arr) {
  if(this->Da != nullptr) {
    delete[] this->Da;
  }
  this->Da = arr;
}

void UHFWfn::setdensb(double *arr) {
  if(this->Db != nullptr) {
    delete[] this->Db;
  }
  this->Db = arr;
}

void UHFWfn::setfocka(double *arr) {
  if(this->Fa != nullptr) {
    delete[] this->Fa;
  }
  this->Fa = arr;
}

void UHFWfn::setfockb(double *arr) {
  if(this->Fb != nullptr) {
    delete[] this->Fb;
  }
  this->Fb = arr;
}

double UHF::energy(const Molecule *molecule,
		      const Wavefunction *wfn_in,
		   UHFWfn *out, std::vector<int> occupieda,
		   std::vector<int> occupiedb) const {

  const SCFWfn *wfn;

  try {
    wfn = static_cast<const SCFWfn *>(wfn_in);
  } catch(std::exception &e) {
    throw(e);
  }

  if(occupieda.size() != 0 && occupieda.size() != (wfn->getelectrons() +
						   wfn->getmultiplicity() - 1)
     / 2) {
    throw(*new std::invalid_argument("Alpha occupation must have the same number of entries as spin-up electrons"));
  }
  if(occupiedb.size() != 0 && occupiedb.size() != (wfn->getelectrons() -
						   wfn->getmultiplicity() + 1)
     / 2) {
    throw(*new std::invalid_argument("Beta occupation must have the same number of entries as spin-up electrons"));
  }

  double enuc = nuclear_repulsion(*molecule);

  const double *S = wfn->getoverlap();
  double *eigreal = new double[wfn->getnorbs()];
  double *eigimag = new double[wfn->getnorbs()];
  double *eigvecs = new double[wfn->getnorbs() * wfn->getnorbs()];
  int orbs = wfn->getnorbs(), elecs = wfn->getelectrons();
  int dim;
  int occa = (wfn->getelectrons() + wfn->getmultiplicity() - 1) / 2;
  int occb = (wfn->getelectrons() - wfn->getmultiplicity() + 1) / 2;
  double e0 = 0, e1 = 0;

  // Set up new matrices.
  double *Da0 = new double[orbs * orbs], *Da1 = new double[orbs * orbs];
  double *Ha = new double[orbs * orbs];
  double *Ca = new double[orbs * orbs];
  double *Fa = new double[orbs * orbs];

  double *Db0 = new double[orbs * orbs], *Db1 = new double[orbs * orbs];
  double *Hb = new double[orbs * orbs];
  double *Cb = new double[orbs * orbs];
  double *Fb = new double[orbs * orbs];

  // Create the hamiltonian and initial Fock matrix.
  const double *T = wfn->getkinetic();
  const double *V = wfn->getpotential();
  const TEIArray *tei = wfn->gettei();

  for(int i = 0; i < orbs * orbs; i++) {
    Ha[i] = T[i] + V[i];
    Fa[i] = T[i] + V[i];
    Hb[i] = T[i] + V[i];
    Fb[i] = T[i] + V[i];
  }

#ifdef NDEBUG
  fprintf(stderr, "Hamiltonian matrix.\n");
  for(int i = 0; i < orbs; i++) {
    fprintf(stderr, "[ ");
    for(int j = 0; j < orbs; j++) {
      fprintf(stderr, "%lf ", H[i * orbs + j]);
    }
    fprintf(stderr, "]\n");
  }
#endif
  
  double *work1 = new double[orbs * orbs],
    *work2 = new double[orbs * orbs];
  int cycles = 0,
    max_cycles = this->getoptions().getintoption("max scf cycles");

  double rms,
    rmsconv = this->getoptions().getdoubleoption("scf rms convergence"),
    enconv = this->getoptions().getdoubleoption("scf energy convergence");

  // Find S^{-1/2}
  double *Shalf = new double[orbs * orbs];
  
  memcpy(work1, S, orbs * orbs * sizeof(double));
  int err = LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'N', 'V', orbs,
			work1, orbs, eigreal, eigimag, nullptr, orbs,
			eigvecs, orbs);

#ifdef NDEBUG
  if(err != 0) {
    fprintf(stderr, "Error computing eigenvalues! Error code %d\n", err);
  }
  fprintf(stderr, "Overlap Eigenvalues.\n");
  for(int i = 0; i < orbs; i++) {
    fprintf(stderr, "%f ", eigreal[i]);
  }
  perror("");
#endif
  
  memset(work1, 0, orbs * orbs * sizeof(double));
  for(int i = 0; i < orbs; i++) {
    if(fabs(eigimag[i]) > 1e-6) {
      fprintf(stderr, "Non-real eigenvalues!\n");
    }
    work1[i * orbs + i] = 1 / sqrt(eigreal[i]);
  }
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, orbs, orbs, orbs,
	      1, eigvecs, orbs, work1, orbs, 0, work2, orbs);
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, orbs, orbs, orbs,
	      1, work2, orbs, eigvecs, orbs, 0, Shalf, orbs);
#ifdef NDEBUG
  fprintf(stderr, "Shalf matrix.\n");
  for(int i = 0; i < orbs; i++) {
    fprintf(stderr, "[ ");
    for(int j = 0; j < orbs; j++) {
      fprintf(stderr, "%lf ", Shalf[i * orbs + j]);
    }
    fprintf(stderr, "]\n");
  }
#endif

  do {
    
    // Find the alpha coefficients and energies.
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, orbs, orbs, orbs,
		1, Fa, orbs, Shalf, orbs, 0, work1, orbs);
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, orbs, orbs, orbs,
		1, Shalf, orbs, work1, orbs, 0, work2, orbs);

#ifdef NDEBUG
    if(cycles == 0) {
      fprintf(stderr, "Fa' matrix.\n");
      for(int i = 0; i < orbs; i++) {
	fprintf(stderr, "[ ");
	for(int j = 0; j < orbs; j++) {
	  fprintf(stderr, "%lf ", work2[i * orbs + j]);
	}
	fprintf(stderr, "]\n");
      }
    }
#endif
    
    LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'N', 'V', orbs,
		  work2, orbs, eigreal, eigimag, nullptr, orbs,
		  eigvecs, orbs);

#ifdef NDEBUG
    if(cycles == 0) {
      fprintf(stderr, "Ca' matrix presort.\n");
      for(int i = 0; i < orbs; i++) {
	fprintf(stderr, "[ ");
	for(int j = 0; j < orbs; j++) {
	  fprintf(stderr, "%lf ", eigvecs[i * orbs + j]);
	}
	fprintf(stderr, "]\n");
      }
      fprintf(stderr, "Energies presort.\n");
      for(int i = 0; i < orbs; i++) {
	fprintf(stderr, "%lf ", eigreal[i]);
      }
      fprintf(stderr, "\n");
    }
#endif
    // Sort the vectors.
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
      double temp = eigreal[i];
      eigreal[i] = eigreal[minj];
      eigreal[minj] = temp;
      
      // Set the vector.
      for(int j = 0; j < orbs; j++) {
	temp = eigvecs[j * orbs + i];
	eigvecs[j * orbs + i] = eigvecs[j * orbs + minj];
	eigvecs[j * orbs + minj] = temp;
      }
    }

    memcpy(work1, eigvecs, orbs * orbs * sizeof(double));
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, orbs, orbs, orbs,
		1, Shalf, orbs, eigvecs, orbs, 0, Ca, orbs);

#ifdef NDEBUG
    if(cycles == 0) {
      fprintf(stderr, "Ca matrix.\n");
      for(int i = 0; i < orbs; i++) {
	fprintf(stderr, "[ ");
	for(int j = 0; j < orbs; j++) {
	  fprintf(stderr, "%lf ", Ca[i * orbs + j]);
	}
	fprintf(stderr, "]\n");
      }
    }
#endif

    // Compute the densities.
    for(int mu = 0; mu < orbs; mu++) {
      for(int nu = 0; nu <= mu; nu++) {
	Da0[mu * orbs + nu] = Da1[mu * orbs + nu];
	Da1[mu * orbs + nu] = 0;
        if(occupieda.size() == 0) {
	  Da1[mu * orbs + nu] += Ca[mu * orbs + m] * Ca[nu * orbs + m];
	} else {
	  Da1[mu * orbs + nu] += Ca[mu * orbs + occupieda[m]] *
	    Ca[nu * orbs + occupieda[m]];
	}
	Da0[nu * orbs + mu] = Da0[mu * orbs + nu];
	Da1[nu * orbs + mu] = Da1[mu * orbs + nu];
      }
    }

#ifdef NDEBUG
    if(cycles == 0) {
      fprintf(stderr, "Initial Alpha Density matrix.\n");
      for(int i = 0; i < orbs; i++) {
	fprintf(stderr, "[ ");
	for(int j = 0; j < orbs; j++) {
	  fprintf(stderr, "%lf ", Da1[i * orbs + j]);
	}
	fprintf(stderr, "]\n");
      }
    }
#endif

    // Find the beta coefficients and energies.
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, orbs, orbs, orbs,
		1, Fb, orbs, Shalf, orbs, 0, work1, orbs);
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, orbs, orbs, orbs,
		1, Shalf, orbs, work1, orbs, 0, work2, orbs);

#ifdef NDEBUG
    if(cycles == 0) {
      fprintf(stderr, "Fb' matrix.\n");
      for(int i = 0; i < orbs; i++) {
	fprintf(stderr, "[ ");
	for(int j = 0; j < orbs; j++) {
	  fprintf(stderr, "%lf ", work2[i * orbs + j]);
	}
	fprintf(stderr, "]\n");
      }
    }
#endif
    
    LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'N', 'V', orbs,
		  work2, orbs, eigreal, eigimag, nullptr, orbs,
		  eigvecs, orbs);

#ifdef NDEBUG
    if(cycles == 0) {
      fprintf(stderr, "Cb' matrix presort.\n");
      for(int i = 0; i < orbs; i++) {
	fprintf(stderr, "[ ");
	for(int j = 0; j < orbs; j++) {
	  fprintf(stderr, "%lf ", eigvecs[i * orbs + j]);
	}
	fprintf(stderr, "]\n");
      }
      fprintf(stderr, "Energies presort.\n");
      for(int i = 0; i < orbs; i++) {
	fprintf(stderr, "%lf ", eigreal[i]);
      }
      fprintf(stderr, "\n");
    }
#endif
    // Sort the vectors.
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
      double temp = eigreal[i];
      eigreal[i] = eigreal[minj];
      eigreal[minj] = temp;
      
      // Set the vector.
      for(int j = 0; j < orbs; j++) {
	temp = eigvecs[j * orbs + i];
	eigvecs[j * orbs + i] = eigvecs[j * orbs + minj];
	eigvecs[j * orbs + minj] = temp;
      }
    }

    memcpy(work1, eigvecs, orbs * orbs * sizeof(double));
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, orbs, orbs, orbs,
		1, Shalf, orbs, eigvecs, orbs, 0, Cb, orbs);

#ifdef NDEBUG
    if(cycles == 0) {
      fprintf(stderr, "Cb matrix.\n");
      for(int i = 0; i < orbs; i++) {
	fprintf(stderr, "[ ");
	for(int j = 0; j < orbs; j++) {
	  fprintf(stderr, "%lf ", Cb[i * orbs + j]);
	}
	fprintf(stderr, "]\n");
      }
    }
#endif

    // Compute the densities.
    for(int mu = 0; mu < orbs; mu++) {
      for(int nu = 0; nu <= mu; nu++) {
	Db0[mu * orbs + nu] = Db1[mu * orbs + nu];
	Db1[mu * orbs + nu] = 0;
        if(occupiedb.size() == 0) {
	  Db1[mu * orbs + nu] += Cb[mu * orbs + m] * Cb[nu * orbs + m];
	} else {
	  Db1[mu * orbs + nu] += Cb[mu * orbs + occupiedb[m]] *
	    Cb[nu * orbs + occupiedb[m]];
	}
	Db0[nu * orbs + mu] = Db0[mu * orbs + nu];
	Db1[nu * orbs + mu] = Db1[mu * orbs + nu];
      }
    }

#ifdef NDEBUG
    if(cycles == 0) {
      fprintf(stderr, "Initial Beta Density matrix.\n");
      for(int i = 0; i < orbs; i++) {
	fprintf(stderr, "[ ");
	for(int j = 0; j < orbs; j++) {
	  fprintf(stderr, "%lf ", Db1[i * orbs + j]);
	}
	fprintf(stderr, "]\n");
      }
    }
#endif

    // Compute the energy.
    e0 = e1;
    e1 = enuc;
    for(int mu = 0; mu < orbs; mu++) {
      for(int nu = 0; nu < orbs; nu++) {
	e1 += Da1[mu * orbs + nu] * (Ha[mu * orbs + nu] + Fa[mu * orbs + nu]);
	e1 += Db1[mu * orbs + nu] * (Hb[mu * orbs + nu] + Fb[mu * orbs + nu]);
      }
    }

    // Compute the new Fock matrix.
    for(int mu = 0; mu < orbs; mu++) {
      for(int nu = 0; nu < orbs; nu++) {
	Fa[mu * orbs + nu] = Ha[mu * orbs + nu];
	Fa[mu * orbs + nu] = Ha[mu * orbs + nu];
	for(int lam = 0; lam < orbs; lam++) {
	  for(int sig = 0; sig < orbs; sig++) {
	    Fa[mu * orbs + nu] += (Da1[lam * orbs + sig] + Db1[lam * orbs + sig]) *
	      tei->at(mu, nu, lam, sig) -
	      Da1[lam * orbs + sig] * tei->at(mu, lam, nu, sig);
	    Fb[mu * orbs + nu] += (Da1[lam * orbs + sig] + Db1[lam * orbs + sig]) *
	      tei->at(mu, nu, lam, sig) -
	      Db1[lam * orbs + sig] * tei->at(mu, lam, nu, sig);
	  }
	}
      }
    }

#ifdef NDEBUG
    if(cycles == 0) {
      fprintf(stderr, "New Fock matrix.\n");
      for(int i = 0; i < orbs; i++) {
	fprintf(stderr, "[ ");
	for(int j = 0; j < orbs; j++) {
	  fprintf(stderr, "%lf ", F[i * orbs + j]);
	}
	fprintf(stderr, "]\n");
      }
    }
#endif

    // Compute the differences.
    rms = 0;
    for(int mu = 0; mu < orbs; mu++) {
      for(int nu = 0; nu < orbs; nu++) {
	rms += (D1[mu * orbs + nu] - D0[mu * orbs + nu]) *
	  (D1[mu * orbs + nu] - D0[mu * orbs + nu]);
      }
    }
    rms /= orbs * orbs;

#ifdef NDEBUG
    if(cycles == 0) {
      fprintf(stderr, "Iteration\tEnergy\t\tDelta\t\tRMS\n");
    }
    if(cycles < 3) {
      fprintf(stderr, "%d\t\t%f\n", cycles, e1);
    } else {
      fprintf(stderr, "%d\t\t%f\t%f\t%f\n", cycles, e1, fabs(e1 - e0), rms);
    }
#endif

    cycles++;
  } while(cycles < max_cycles &&
	  (cycles < 3 || rms > rmsconv || fabs(e1 - e0) > enconv));

#ifdef NDEBUG
  fprintf(stderr, "Final Fock matrix.\n");
  for(int i = 0; i < orbs; i++) {
    fprintf(stderr, "[ ");
    for(int j = 0; j < orbs; j++) {
      fprintf(stderr, "%lf ", F[i * orbs + j]);
    }
    fprintf(stderr, "]\n");
  }
#endif

  // Set the wavefunction return.
  if(out != nullptr) {
    out->setcoef(C);
    out->setdens(D1);
    out->setfock(F);
    out->setenergies(eigreal);
  } else {
    delete[] C;
    delete[] D1;
    delete[] F;
    delete[] eigreal;
  }
  delete[] eigimag;
  delete[] eigvecs;
  delete[] D0;
  delete[] H;
  delete[] work1;
  delete[] work2;
  delete[] Shalf;

  if(cycles >= max_cycles && (rms > rmsconv || fabs(e1 - e0) > enconv)) {
    throw(*new std::runtime_error("SCF energy did not converge!"));
  }

  return e1;
}
  
double UHF::energy(const Molecule *mol,
		      const Wavefunction *wfn_in) const {
  return this->energy(mol, wfn_in, nullptr);
}
