#ifndef INTEGRAL_FACTORY_CPP
#define INTEGRAL_FACTORY_CPP


#include "integrals.hpp"
#include "../util/molecule.hpp"
#include "../util/tei_array.hpp"
#include "../util/atom.hpp"
#include <exception>
#include <stdexcept>
#include <vector>
#include <thread>
#include <cmath>
#include "../opts/options.hpp"

template<typename Ints>
compchem::IntegralFactory<Ints>::~IntegralFactory() {
  delete &(this->opts);
}

template<typename Ints>
void compchem::IntegralFactory<Ints>::s_routine(const std::vector<const compchem::GaussianOrbital *> *orbs,
		      const std::vector<std::array<double, 3> *> *centers,
		      double *out,
		      int dim, int thread_num, int threads,
		      compchem::OptionList &opts) {
  Ints method = Ints(opts);

  int max_size = (dim * (dim + 1)) / 2;

  for(int i = (int) (thread_num * (double) max_size / threads);
      i < ((thread_num == threads - 1)? max_size:
	   (int) ((thread_num + 1) * (double) max_size / threads));
      i++) {
    // Calculate the indices.
    int mu = (int) std::floor((std::sqrt(1 + 8 * i) - 1) / 2.0);
    int nu = i - (mu * (mu + 1)) / 2;

    // Compute.
    double res = method.overlap(orbs->at(mu),
				orbs->at(nu),
				*centers->at(mu),
				*centers->at(nu));
    out[mu * dim + nu] = res;
    out[nu * dim + mu] = method.overlap(orbs->at(nu),
				orbs->at(mu),
				*centers->at(nu),
				*centers->at(mu));
  }
}

template<typename Ints>
void compchem::IntegralFactory<Ints>::t_routine(const std::vector<const compchem::GaussianOrbital *> *orbs,
		      const std::vector<std::array<double, 3> *> *centers,
		      double *out,
		      int dim, int thread_num, int threads,
		      compchem::OptionList &opts) {
  Ints method = Ints(opts);

  int max_size = (dim * (dim + 1)) / 2;

  for(int i = (int) (thread_num * (double) max_size / threads);
      i < ((thread_num == threads - 1)? max_size:
	   (int) ((thread_num + 1) * (double) max_size / threads));
      i++) {
    // Calculate the indices.
    int mu = (int) std::floor((std::sqrt(1 + 8 * i) - 1) / 2.0);
    int nu = i - (mu * (mu + 1)) / 2;

    // Compute.
    double res = method.laplacian(orbs->at(mu),
				orbs->at(nu),
				*centers->at(mu),
				*centers->at(nu));
    out[mu * dim + nu] = res;
    out[nu * dim + mu] = res;
  }
}

template<typename Ints>
void compchem::IntegralFactory<Ints>::v_routine(const std::vector<const compchem::GaussianOrbital *> *orbs,
						const std::vector<std::array<double, 3> *> *centers,
						double *out,
						int dim, int thread_num, int threads,
						compchem::OptionList &opts,
						const compchem::Molecule *mol) {
  Ints method = Ints(opts);

  int max_size = (dim * (dim + 1)) / 2;

  for(int i = (int) (thread_num * (double) max_size / threads);
      i < ((thread_num == threads - 1)? max_size:
	   (int) ((thread_num + 1) * (double) max_size / threads));
      i++) {
    // Calculate the indices.
    int mu = (int) std::floor((std::sqrt(1 + 8 * i) - 1) / 2.0);
    int nu = i - (mu * (mu + 1)) / 2;
    double sum = 0;
    for(int j = 0; j < mol->getsize(); j++) {
      // Compute.
      sum += method.coulomb(orbs->at(mu),
			    orbs->at(nu),
			    *centers->at(mu),
			    *centers->at(nu),
			    mol->getatom(j));
    }
    out[mu * dim + nu] = sum;
    out[nu * dim + mu] = sum;
  }
}

template<typename Ints>
void compchem::IntegralFactory<Ints>::tei_routine(const std::vector<const compchem::GaussianOrbital *> *orbs,
		      const std::vector<std::array<double, 3> *> *centers,
		      TEIArray *out,
		      int dim, int thread_num, int threads,
		      compchem::OptionList &opts) {
  Ints method = Ints(opts);

  int max_size = out->getsize();

  for(int i = (int) (thread_num * (double) max_size / threads);
      i < ((thread_num == threads - 1)? max_size:
	   (int) ((thread_num + 1) * (double) max_size / threads));
      i++) {
    int mu, nu, lam, sig;
    out->indextoquad(i, &mu, &nu, &lam, &sig);
    out->at_direct(i) = method.exchange(orbs->at(mu),
					orbs->at(nu),
					orbs->at(lam),
					orbs->at(sig),
					*centers->at(mu),
					*centers->at(nu),
					*centers->at(lam),
					*centers->at(sig));
  }
}
  
  

template<typename Ints>
void compchem::IntegralFactory<Ints>::Smatrix(const compchem::Molecule *mol, double *out, int *dim) {
  if(out == nullptr || dim == nullptr) {
    throw new std::invalid_argument("The output values to the integral factory cn not be null.");
  }
  if(mol == nullptr) {
    throw new std::invalid_argument("The molecule input to the integral factory can not be null.");
  }
  
  int orbs = 0;
  std::vector<const compchem::GaussianOrbital *> *orbitals =
    new std::vector<const compchem::GaussianOrbital *>();
  std::vector<std::array<double, 3> *> *centers =
    new std::vector<std::array<double, 3> *>();
  for(int i = 0; i < mol->getsize(); i++) {
    orbs += mol->getatom(i).getnorbitals();
    for(int j = 0; j < mol->getatom(i).getnorbitals(); j++) {
      orbitals->push_back(static_cast<const compchem::GaussianOrbital *>
			 (&mol->getatom(i).getorbital(j)));
      centers->push_back(new std::array<double, 3>(mol->getatom(i).getpos()));
    }
  }
  *dim = orbs;

  std::vector<std::thread> threads;

  int nthreads = this->opts.getintoption("threads");
  if(nthreads > (*dim * (*dim + 1)) / 2) {
    nthreads = (*dim * (*dim + 1)) / 2;
  }

  for(int i = 0; i < nthreads - 1; i++) {
    threads.push_back(std::thread(IntegralFactory<Ints>::s_routine, orbitals,
				  centers,
				  out, *dim, i,
				  nthreads,
				  std::ref(this->opts)));
  }

  IntegralFactory<Ints>::s_routine(orbitals, centers, out, *dim,
				   nthreads - 1,
				   nthreads,
				   this->opts);
  
  for(int i = 0; i < threads.size(); i++) {
    threads[i].join();
  }

  for(int i = 0; i < centers->size(); i++) {
    delete centers->at(i);
  }
  orbitals->clear();
  centers->clear();

  delete orbitals;
  delete centers;
  
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
  std::vector<const compchem::GaussianOrbital *> *orbitals =
    new std::vector<const compchem::GaussianOrbital *>();
  std::vector<std::array<double, 3> *> *centers =
    new std::vector<std::array<double, 3> *>();
  for(int i = 0; i < mol->getsize(); i++) {
    orbs += mol->getatom(i).getnorbitals();
    for(int j = 0; j < mol->getatom(i).getnorbitals(); j++) {
      orbitals->push_back(static_cast<const compchem::GaussianOrbital *>
			 (&mol->getatom(i).getorbital(j)));
      centers->push_back(new std::array<double, 3>(mol->getatom(i).getpos()));
    }
  }
  *dim = orbs;

  std::vector<std::thread> threads;

  int nthreads = this->opts.getintoption("threads");
  if(nthreads > (*dim * (*dim + 1)) / 2) {
    nthreads = (*dim * (*dim + 1)) / 2;
  }

  for(int i = 0; i < nthreads - 1; i++) {
    threads.push_back(std::thread(IntegralFactory<Ints>::t_routine, orbitals,
				  centers,
				  out, *dim, i,
				  nthreads,
				  std::ref(this->opts)));
  }

  IntegralFactory<Ints>::t_routine(orbitals, centers, out, *dim,
				   nthreads - 1,
				   nthreads,
				   this->opts);
  
  for(int i = 0; i < threads.size(); i++) {
    threads[i].join();
  }

  for(int i = 0; i < centers->size(); i++) {
    delete centers->at(i);
  }
  orbitals->clear();
  centers->clear();

  delete orbitals;
  delete centers;
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
  std::vector<const compchem::GaussianOrbital *> *orbitals =
    new std::vector<const compchem::GaussianOrbital *>();
  std::vector<std::array<double, 3> *> *centers =
    new std::vector<std::array<double, 3> *>();
  for(int i = 0; i < mol->getsize(); i++) {
    orbs += mol->getatom(i).getnorbitals();
    for(int j = 0; j < mol->getatom(i).getnorbitals(); j++) {
      orbitals->push_back(static_cast<const compchem::GaussianOrbital *>
			 (&mol->getatom(i).getorbital(j)));
      centers->push_back(new std::array<double, 3>(mol->getatom(i).getpos()));
    }
  }
  *dim = orbs;

  for(int i = 0; i < orbs * orbs; i++) {
    out[i] = 0;
  }

  std::vector<std::thread> threads;
  
  int nthreads = this->opts.getintoption("threads");
  if(nthreads > (*dim * (*dim + 1)) / 2) {
    nthreads = (*dim * (*dim + 1)) / 2;
  }

  for(int i = 0; i < nthreads - 1; i++) {
    threads.push_back(std::thread(IntegralFactory<Ints>::v_routine, orbitals,
				  centers,
				  out, *dim, i,
				  nthreads,
				  std::ref(this->opts), mol));
  }

  IntegralFactory<Ints>::v_routine(orbitals, centers, out, *dim,
				   nthreads - 1,
				   nthreads,
				   this->opts, mol);

  for(int i = 0; i < threads.size(); i++) {
    threads[i].join();
  }

  for(int i = 0; i < centers->size(); i++) {
    delete centers->at(i);
  }
  orbitals->clear();
  centers->clear();

  delete orbitals;
  delete centers;
}

template<typename Ints>
compchem::TEIArray *compchem::IntegralFactory<Ints>::TEIints(const compchem::Molecule *mol) {
  if(mol == nullptr) {
    throw new std::invalid_argument("The molecule input to the integral factory can not be null.");
  }
  
  int orbs = 0;
  std::vector<const compchem::GaussianOrbital *> *orbitals =
    new std::vector<const compchem::GaussianOrbital *>();
  std::vector<std::array<double, 3> *> *centers =
    new std::vector<std::array<double, 3> *>();
  for(int i = 0; i < mol->getsize(); i++) {
    orbs += mol->getatom(i).getnorbitals();
    for(int j = 0; j < mol->getatom(i).getnorbitals(); j++) {
      orbitals->push_back(static_cast<const compchem::GaussianOrbital *>
			 (&mol->getatom(i).getorbital(j)));
      centers->push_back(new std::array<double, 3>(mol->getatom(i).getpos()));
    }
  }

  TEIArray *out = new TEIArray(orbs);

  std::vector<std::thread> threads;

  int nthreads = this->opts.getintoption("threads");
  if(nthreads > out->getsize()) {
    nthreads = out->getsize();
  }

  for(int i = 0; i < nthreads - 1; i++) {
    threads.push_back(std::thread(IntegralFactory<Ints>::tei_routine, orbitals,
				  centers,
				  out, orbs, i,
				  nthreads,
				  std::ref(this->opts)));
  }

  IntegralFactory<Ints>::tei_routine(orbitals, centers, out, orbs,
				   nthreads - 1,
				   nthreads,
				   this->opts);

  for(int i = 0; i < threads.size(); i++) {
    threads[i].join();
  }

  for(int i = 0; i < centers->size(); i++) {
    delete centers->at(i);
  }
  orbitals->clear();
  centers->clear();

  delete orbitals;
  delete centers;

  return out;
}
  

#endif
