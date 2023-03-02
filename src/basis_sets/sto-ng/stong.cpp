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
  // Beautiful identity that I found.
  double arit = (alpha1 + alpha2) / 2,
    geom = std::sqrt(alpha1 * alpha2);
  return std::pow(geom / arit, l + 1.5);
}

double compchem::STO_nGComputer::gderiva1(int l, double alpha1,
					  double alpha2) const {
  return this->gfunc(l, alpha1, alpha2) *
    ((2 * l + 3) / (4 * alpha1) - (l + 1.5) / (alpha1 + alpha2));
}

double compchem::STO_nGComputer::sfunc(int N, int l, double alpha,
				      double zeta) const {
  return 2 * std::pow(8 * alpha * alpha * alpha / M_PI, 0.25) *
    std::pow(2 * zeta, N) *
    std::sqrt(std::pow(8 * alpha, l) * zeta *
	      std::tgamma(l + 2) / std::tgamma(2 * l + 3) /
	      std::tgamma(2 * N + 1)) *
    std::pow(alpha, -(N + l + 2.0) / 2.0) *
    (std::tgamma((N + l + 2.0) / 2.0) *
     hyper1f1((N + l + 2.0) / 2.0, 0.5, zeta * zeta / (4 * alpha)) -
     zeta / std::sqrt(alpha) * std::tgamma((N + l + 3.0) / 2.0) *
     hyper1f1((N + l + 3.0) / 2.0, 1.5, zeta * zeta / (4 * alpha)));
}

double compchem::STO_nGComputer::sderiv(int N, int l, double alpha,
				      double zeta) const {
  // Until I can find the analytic derivative, use the numeric derivative.
  return (this->sfunc(N, l, alpha + 0.00001 * alpha, zeta) -
	  this->sfunc(N, l, alpha, zeta)) / (0.00001 * alpha);

  
  /*
  return -(2 * N + 1.0) / (4 * alpha) * this->sfunc(N, l, alpha, zeta) +
    2 * std::pow(8 * alpha * alpha * alpha / M_PI, 0.25) *
    std::pow(2 * zeta, N) *
    std::sqrt(zeta * std::pow(8 * alpha, l) * std::tgamma(l + 2) /
	      std::tgamma(2 * N + 1) / std::tgamma(2 * l + 3)) *
    std::pow(alpha, -(N + l + 7.0) / 2.0) * zeta / 12.0 *
    (zeta * zeta * (N + l + 3) * std::tgamma((N + l + 3.0) / 2.0) *
     hyper1f1((N + l + 5.0) / 2.0, 2.5, zeta * zeta / (4 * alpha)) -
     3 * std::sqrt(alpha) * zeta * (N + l + 2) *
     std::tgamma((N + l + 3.0) / 2.0) *
     hyper1f1((N + l + 4.0) / 2.0, 1.5, zeta * zeta / (4 * alpha)) +
     6 * alpha * std::tgamma((N + l + 3.0) / 2.0) *
     hyper1f1((N + l + 3.0) / 2.0, 1.5, zeta * zeta / (4 * alpha)));
  */
}


std::vector<double> compchem::STO_nGComputer::gradient(
				       const std::vector<double> &coefs,
				       const std::vector<double> &alphas,
				       double zeff, int N, int l,
				       double lambda) const {
  
  std::vector<double> out = std::vector<double>(2 * this->getnGauss() + 1);

  for(int i = 0; i < out.size(); i++) {
    out[i] = 0;
  }

  // Compute the coefficient gradient.
  for(int i = 0; i < this->getnGauss(); i++) {
    double sum1 = 2 * coefs[i];
    for(int j = 0; j < this->getnGauss(); j++) {
      if(i == j) {
	continue;
      }
      sum1 += coefs[j] * this->gfunc(l, alphas[i], alphas[j]);
    }
    out[i] = -this->sfunc(N, l, alphas[i], zeff) - lambda * sum1;
  }
  // Compute the alpha gradient.
  for(int i = 0; i < this->getnGauss(); i++) {
    double sum1 = 0;
    for(int j = 0; j < this->getnGauss(); j++) {
      if(i == j) {
	// Not needed, but it should avoid computing trivial values.
	continue;
      }
      sum1 += coefs[i] * coefs[j] * this->gderiva1(l, alphas[i], alphas[j]);
    }
    out[i + this->getnGauss()] = -coefs[i] *
      this->sderiv(N, l, alphas[i], zeff) - lambda * sum1;
  }

  // Compute the contribution of the Lagrange multiplier.
  out[2 * this->getnGauss()] = 1;
  for(int i = 0; i < this->getnGauss(); i++) {
    out[2 * this->getnGauss()] -= coefs[i] * coefs[i];
    for(int j = 0; j < i; j++) {
      out[2 * this->getnGauss()] -= 2 * coefs[i] * coefs[j] * this->gfunc(l, alphas[i], alphas[j]);
    }
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
    sum -= lambda * coefs[i] * coefs[i];
  }
  sum += 1;   // Slater

  // Compute the Gaussian-Gaussian term.
  for(int i = 0; i < this->getnGauss(); i++) {
    for(int j = 0; j < i; j++) {
      sum -= lambda * 2 * coefs[i] * coefs[j] *
	this->gfunc(l, alphas[i], alphas[j]);
    }
  }

  // Compute the Gaussian-Slater term.
  for(int i = 0; i < this->getnGauss(); i++) {
    sum += this->sfunc(N, l, alphas[i], zeff) * coefs[i];
  }

  return sum;
}
    
      

#define CONVERGENCE 1e-8
#define STEP 1e-2
#define MAX_CYCLES 1000
void compchem::STO_nGComputer::gradientdescent(
		       std::vector<double> &coefs,
		       std::vector<double> &alphas,
		       double zeff, int N, int l) const {
  
  std::vector<double> grad1 = std::vector<double>(2 * this->getnGauss() + 1),
    grad2 = std::vector<double>(2 * this->getnGauss() + 1),
    coefs2 = std::vector<double>(this->getnGauss()),
    alphas2 = std::vector<double>(this->getnGauss());
  double step, lambda1 = 1, lambda2;
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
    std::fprintf(stderr, "], Lambda: %lf.\n", lambda1);
    std::fprintf(stderr, "Loss: %lf.\n\n", this->loss(coefs, alphas, zeff, N, l, lambda1));

    alphas2 = alphas;
    coefs2 = coefs;
    lambda2 = lambda1;
    
    // Descend
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

std::vector<compchem::GaussianOrbital> *compchem::STO_nGComputer::compute(
        int Z, const std::vector<double> &guess_c,
	const std::vector<double> &guess_a) const {
  
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
      if(i < guess_c.size()) {
	coefs->at(i) = guess_c[i];
      } else if(i == 0) {
	coefs->at(i) = 1;
      } else {
	coefs->at(i) = 0;
      }
      if(i < guess_a.size()) {
	alphas->at(i) = guess_a[i];
      } else {
	alphas->at(i) = 8 * std::sqrt(zeff) / (i + 1);
      }
    }

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
      if(i == 0) {
	coefs->at(i) = 1;
      } else {
	coefs->at(i) = 0;
      }
      alphas->at(i) = 8 * std::sqrt(zeff) / (i + 1);
    }

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
