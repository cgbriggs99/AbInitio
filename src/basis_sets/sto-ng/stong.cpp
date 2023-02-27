/*
 * stong.cpp
 *
 *  Created on: Feb 7, 2023
 *      Author: connor
 */

#include "stong.hpp"
#include "../../util/electron_configs.hpp"
#include <array>
#include <cmath>
#include <extramath.h>
#include <vector>
#include "../../util/vector_ops.hpp"
#include <cstdlib>
#include <cstdio>

int compchem::STO_nGComputer::getnGauss() const {
	return this->n;
}

double compchem::STO_nGComputer::gfunc(int l, double alpha1,
				       double alpha2) const {
  return std::pow(4 * alpha1 * alpha2 / M_PI / M_PI, 0.75) * std::tgamma(l + 1.5) /
    std::pow(alpha1 + alpha2, l + 1.5) *
    M_PI * std::pow(64 * alpha1 * alpha2, l / 2.0) *
    std::exp(std::lgamma(l + 2) - std::lgamma(2 * l + 3)) * 4;
}

double compchem::STO_nGComputer::gderiva1(int l, double alpha1,
					  double alpha2) const {
  return this->gfunc(l, alpha1, alpha2) *
    (l / 2.0 / alpha1 + 3.0 / (4.0 * alpha1) - (l + 1.5) / (alpha1 + alpha2)) * 4;
}

double compchem::STO_nGComputer::sfunc(int N, int l, double alpha,
				      double zeta) const {
  return std::pow(2 * alpha / M_PI, 0.75) * std::pow(2 * zeta, N) *
    std::sqrt(2 * zeta / std::tgamma(2 * N + 1)) *
    std::pow(alpha, -(N + l + 2.0) / 2.0) * (std::tgamma((N + l + 2.0) / 2.0) *
	   hyper1f1((N + l + 2.0) / 2.0, 0.5, zeta * zeta / 4 / alpha) -
	   zeta / std::sqrt(alpha) * std::tgamma((N + l + 3.0) / 2.0) *
	   hyper1f1((N + l + 3.0) / 2.0, 1.5, zeta * zeta / 4 / alpha)) * 2 *
    std::sqrt(M_PI * std::pow(8 * alpha, l) * std::exp(std::lgamma(l + 2) -
						       std::lgamma(2 * l + 3)));
}

double compchem::STO_nGComputer::sderiv(int N, int l, double alpha,
				      double zeta) const {
  return 2 * std::pow(2 * zeta, N) * std::sqrt(2 * zeta /
						 std::tgamma(2 * N + 1)) *
    std::pow(alpha, -(2 * N + 2 * l + 11.0) / 4.0) / (6 * std::pow(2, 0.25) *
						      std::pow(M_PI, 0.75)) *
    (zeta * std::tgamma((N + l + 3.0) / 2.0) *
     (zeta * zeta * (N + l + 3) *
      hyper1f1((N + l + 5.0) / 2.0, 2.5, zeta * zeta / 4 / alpha) +
      3 * alpha * (2 * N + 2 * l + 3) *
      hyper1f1((N + l + 3.0) / 2.0, 1.5, zeta * zeta / 4 / alpha)) -
     3 * std::sqrt(alpha) * std::tgamma((N + l + 2.0) / 2.0) *
     (zeta * zeta * (N + l + 2) *
      hyper1f1((N + l + 4.0) / 2.0, 1.5, zeta * zeta / 4 / alpha) +
      alpha * (2 * N + 2 * l + 1) *
      hyper1f1((N + l + 2.0) / 2.0, 0.5, zeta * zeta / 4 / alpha))) +
    this->sfunc(N, l, alpha, zeta) / alpha;
}


std::vector<double> compchem::STO_nGComputer::gradient(
				       const std::vector<double> &coefs,
				       const std::vector<double> &alphas,
				       double zeff, int N, int l,
				       double lambda) const {
  
  std::vector<double> out = std::vector<double>(2 * this->getnGauss() + 1);

  // Compute the coefficient gradient.
  for(int i = 0; i < this->getnGauss(); i++) {
    double sum1 = 2 * coefs[i];
    for(int j = 0; j < i; j++) {
      sum1 += 2 * coefs[j] * this->gfunc(l, alphas[i], alphas[j]);
    }
    out[i] = (1 - lambda) * sum1 - this->sfunc(N, l, alphas[i], zeff);
  }
  // Compute the alpha gradient.
  for(int i = 0; i < this->getnGauss(); i++) {
    double sum1 = 0;
    for(int j = 0; j < i; j++) {
      sum1 += 2 * coefs[i] * coefs[j] * this->gderiva1(l, alphas[i], alphas[j]);
    }
    out[i + this->getnGauss()] += sum1 * (1 - lambda) -
      this->sderiv(N, l, alphas[i], zeff);
  }

  // Compute the contribution of the Lagrange multiplier.
  out[2 * this->getnGauss()] = 1;
  for(int i = 0; i < this->getnGauss(); i++) {
    out[2 * this->getnGauss()] -= coefs[i] * coefs[i];
    double sum2 = 0;
    for(int j = 0; j < i; j++) {
      sum2 -= 2 * coefs[i] * coefs[j] * this->gfunc(l, alphas[i], alphas[j]);
    }
    out[2 * this->getnGauss()] -= 2 * coefs[i] * sum2;
  }
    
  return out;
  
}

std::vector<double> compchem::STO_nGComputer::gradientptype(
			       const std::vector<double> &coefs,
			       const std::vector<double> &alphas,
			       double zeff, int N, int l, double lambda) const {
  std::vector<double> out = std::vector<double>(this->getnGauss() + 1);

  // Compute the coefficient gradient.
  for(int i = 0; i < this->getnGauss(); i++) {
    double sum1 = 2 * coefs[i];
    for(int j = 0; j < i; j++) {
      sum1 += 2 * coefs[j] * this->gfunc(l, alphas[i], alphas[j]);
    }
    out[i] = (1 - lambda) * sum1 - this->sfunc(N, l, alphas[i], zeff);
  }

  // Compute the contribution of the Lagrange multiplier.
  out[this->getnGauss()] = 1;
  for(int i = 0; i < this->getnGauss(); i++) {
    out[this->getnGauss()] -= coefs[i] * coefs[i];
    double sum2 = 0;
    for(int j = 0; j < i; j++) {
      sum2 -= 2 * coefs[i] * coefs[j] * this->gfunc(l, alphas[i], alphas[j]);
    }
    out[this->getnGauss()] -= 2 * coefs[i] * sum2;
  }

  return out;
}

  

double compchem::STO_nGComputer::loss(const std::vector<double> &coefs,
			       const std::vector<double> &alphas,
			       double zeff, int N, int l,
			       double lambda) const {
  double sum = lambda;    // Initialize with the Lagrangian multiplier.

  // Compute the orthogonal terms.
  for(int i = 0; i < this->getnGauss(); i++) {    // Gaussian
    sum += coefs[i] * coefs[i] * (-lambda + 1);
  }
  sum += 1;   // Slater

  // Compute the Gaussian-Gaussian term.
  for(int i = 0; i < this->getnGauss(); i++) {
    double sum2 = 0;
    for(int j = 0; j < this->getnGauss(); j++) {
      if(i == j) {
	continue;
      }
      sum2 += coefs[j] * this->gfunc(l, alphas[i], alphas[j]);
    }
    sum += 2 * coefs[i] * sum2 * (-lambda + 1);
  }

  // Compute the Gaussian-Slater term.
  for(int i = 0; i < this->getnGauss(); i++) {
    sum -= this->sfunc(N, l, alphas[i], zeff) * coefs[i];
  }

  return sum;
}
    
      

#define CONVERGENCE 1e-8
#define STEP 1e-3
#define MAX_CYCLES 100
void compchem::STO_nGComputer::gradientdescent(
		       std::vector<double> &coefs,
		       std::vector<double> &alphas,
		       double zeff, int N, int l) const {
  
  std::vector<double> grad1 = std::vector<double>(2 * this->getnGauss() + 1),
    grad2 = std::vector<double>(2 * this->getnGauss() + 1),
    coefs2 = std::vector<double>(this->getnGauss()),
    alphas2 = std::vector<double>(this->getnGauss());
  double step, lambda1 = 2, lambda2;
  int cycle = 0;

  do {
    grad2 = grad1;
    
    grad1 = this->gradient(coefs, alphas, zeff, N, l, lambda1);

    double num = 0, den = 0;
    if(cycle >= 3) {
      for(int i = 0; i < this->getnGauss(); i++) {
	num += (coefs[i] - coefs2[i]) *
	  (grad1[i] - grad2[i]) + (alphas[i] - alphas2[i]) *
	  (grad1[i + this->getnGauss()] - grad2[i + this->getnGauss()]);
	den += (grad1[i] - grad2[i]) * (grad1[i] - grad2[i]) +
	  (grad1[i + this->getnGauss()] - grad2[i + this->getnGauss()]) *
	  (grad1[i + this->getnGauss()] - grad2[i + this->getnGauss()]);
      }
      // Don't forget about lambda.
      num += (lambda1 - lambda2) *
	(grad1[2 * this->getnGauss()] - grad2[2 * this->getnGauss()]);
      den += (grad1[2 * this->getnGauss()] - grad2[2 * this->getnGauss()]) *
	(grad1[2 * this->getnGauss()] - grad2[2 * this->getnGauss()]);
    }
      
    // Compute the step. Have a default value if converged.
    if(cycle < 3 || num == 0 || den == 0) {
      step = STEP;
    } else {
      step = fabs(num) / fabs(den);
    }
    if(!isfinite(step)) {
      step = STEP;
    }
    std::fprintf(stderr, "%d: Taking step of size %lf in direction [ ", cycle, step);
    for(int i = 0; i < 2 * this->getnGauss() + 1; i++) {
      std::fprintf(stderr, "%lf ", grad1[i]);
    }
    std::fprintf(stderr, "]\nFrom Coefficients: [ ");
    for(int i = 0; i < this->getnGauss(); i++) {
      std::fprintf(stderr, "%lf ", coefs[i]);
    }
    std::fprintf(stderr, "], Alphas: [ ");
    for(int i = 0; i < this->getnGauss(); i++) {
      std::fprintf(stderr, "%lf ", alphas[i]);
    }
    std::fprintf(stderr, "], Lambda: %lf.\n\n", lambda1);

    alphas2 = alphas;
    coefs2 = coefs;
    lambda2 = lambda1;
    
    // Descend.
    for(int i = 0; i < this->getnGauss(); i++) {
      coefs[i] -= step * grad1[i];
      alphas[i] -= step * grad1[i + this->getnGauss()];
    }
    lambda1 -= step * grad1[2 * this->getnGauss()];
    cycle++;
    // Check for convergence, but don't exit too early.
  } while(cycle < 10 || (cycle < MAX_CYCLES &&
			 vecmag(grad1) / this->getnGauss() > CONVERGENCE &&
			 fabs(this->loss(coefs, alphas,
					 zeff, N, l, lambda1)) > CONVERGENCE));

  std::fprintf(stderr, "Descended in %d cycles.\n", cycle);
  std::fprintf(stderr, "Lagrangian multiplier: %lf\n", lambda1);
  std::fprintf(stderr, "Loss function: %lf.\n", this->loss(coefs, alphas,
							   zeff, N, l, lambda1));
  std::fprintf(stderr, "Gradient: [ ");
  for(int i = 0; i < 2 * this->getnGauss() + 1; i++) {
    std::fprintf(stderr, "%lf ", grad1[i]);
  }
  std::perror("]");
    
  return;
}

void compchem::STO_nGComputer::gradientdescentptype(
		       std::vector<double> &coefs,
		       const std::vector<double> &alphas,
		       double zeff, int N, int l) const {
  
  std::vector<double> grad1 = std::vector<double>(this->getnGauss() + 1),
    grad2 = std::vector<double>(this->getnGauss() + 1),
    coefs2 = std::vector<double>(this->getnGauss());
  double step, lambda1 = 1, lambda2;
  int cycle = 0;

  do {
    grad2 = grad1;
    
    grad1 = this->gradientptype(coefs, alphas, zeff, N, l, lambda1);

    double num = 0, den = 0;
    if(cycle >= 3) {
      for(int i = 0; i < this->getnGauss(); i++) {
	num += (coefs[i] - coefs2[i]) *
	  (grad1[i] - grad2[i]);
	den += (grad1[i] - grad2[i]) * (grad1[i] - grad2[i]);
      }
      // Don't forget about lambda.
      num += (lambda1 - lambda2) *
	(grad1[this->getnGauss()] - grad2[this->getnGauss()]);
      den += (grad1[this->getnGauss()] - grad2[this->getnGauss()]) *
	(grad1[this->getnGauss()] - grad2[this->getnGauss()]);
    }
      
    // Compute the step. Have a default value if converged.
    if(num == 0 || den == 0 || cycle < 3) {
      step = STEP;
    } else {
      step = fabs(num) / den;
    }
    
    coefs2 = coefs;
    lambda2 = lambda1;
    
    // Descend.
    for(int i = 0; i < this->getnGauss(); i++) {
      coefs[i] -= step * grad1[i];
    }
    lambda1 -= step * grad1[this->getnGauss()];
    cycle++;
    // Check for convergence, but don't exit too early.
  } while(cycle < 10 || (vecmag(grad1) / this->getnGauss() > CONVERGENCE &&
	  rmsnorm(grad1, grad2) > CONVERGENCE &&
		    this->loss(coefs, alphas,
			       zeff, N, l, lambda1) > CONVERGENCE));
  return;
}

std::vector<compchem::GaussianOrbital> *compchem::STO_nGComputer::compute(int Z) const {
  static std::array<int, 19> lvals = std::array<int, 19>({0, 0, 1, 0, 1, 0, 2, 1, 0, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1});
  static std::array<int, 19> nvals = std::array<int, 19>({1, 2, 2, 3, 3, 4, 3, 4, 5, 4, 5, 6, 4, 5, 6, 7, 5, 6, 7});
  
  std::vector<GaussianOrbital> *out = new std::vector<GaussianOrbital>();

  std::vector<double> *last_stype_alpha = nullptr;
  GSConfig *config = compchem::getconfig(Z);
  
  for(int i = 0; i < config->getShells(); i++) {
    // Compute the Slater parameter.
    double zeff = Z - slater_rule(nvals[i], lvals[i], *config);

    // Compute the coefficients and scales.
    std::vector<double> *coefs = new std::vector<double>(this->getnGauss()),
      *alphas = new std::vector<double>(this->getnGauss());

    // Set initial guess.
    for(int i = 0; i < this->getnGauss(); i++) {
      coefs->at(i) = 0;
      alphas->at(i) = std::cos(i) * std::cos(i);
    }
    coefs->at(0) = 1;

    // s-type orbital, so store a pointer to it.
    if(lvals[i] == 0) {
      this->gradientdescent(*coefs, *alphas, zeff, nvals[i], lvals[i]);
      if(last_stype_alpha != nullptr) {
	delete last_stype_alpha;
	last_stype_alpha = nullptr;
      }
      last_stype_alpha = new std::vector<double>(*alphas);
    } else if(lvals[i] == 1) {    // p-type orbital, so use s-type alphas.
      this->gradientdescentptype(*coefs, *last_stype_alpha, zeff, nvals[i],
				 lvals[i]);
    } else {
      this->gradientdescent(*coefs, *alphas, zeff, nvals[i], lvals[i]);
    }

    // Add all the orbitals.
    for(int ml = -lvals[i]; ml <= lvals[i]; ml++) {
      if(lvals[i] == 1) {
	out->push_back(GaussianOrbital(lvals[i], ml, *coefs,
					   *last_stype_alpha));
      } else {
	out->push_back(GaussianOrbital(lvals[i], ml, *coefs, *alphas));
      }
    }
    delete coefs;
    delete alphas;
  }
  delete config;
  if(last_stype_alpha != nullptr) {
    delete last_stype_alpha;
  }

  return out;
  
}
