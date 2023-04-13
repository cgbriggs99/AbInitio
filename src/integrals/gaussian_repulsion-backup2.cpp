#include <extramath.h>
#include "../basis_sets/basis_set.hpp"
#include "../util/atom.hpp"
#include <array>
#include <cmath>
#include <vector>
#include <deque>
#include <map>
#include "../util/polynomial.hpp"
#include "integrals.hpp"

using namespace compchem;

static double compute_abcd(const std::array<int, 12> &index,
			   std::map<std::array<int, 12>, double> &abcd,
			   std::map<std::array<int, 12>, double> &apcq,
			   std::map<std::array<int, 4>, double> &rm,
			   const double *f0,
			   double alpha, double beta,
			   double gamma, double delta,
			   double Rx, double Ry, double Rz,
			   const std::array<double, 3> &c1,
			   const std::array<double, 3> &c2,
			   const std::array<double, 3> &c3,
			   const std::array<double, 3> &c4);

static double compute_apcq(const std::array<int, 12> &index,
			   std::map<std::array<int, 12>, double> &apcq,
			   std::map<std::array<int, 4>, double> &rm,
			   const double *f0,
			   double alpha, double beta,
			   double gamma, double delta,
			   double Rx, double Ry, double Rz,
			   const std::array<double, 3> &c1,
			   const std::array<double, 3> &c2,
			   const std::array<double, 3> &c3,
			   const std::array<double, 3> &c4);

static double compute_rm(const std::array<int, 4> &index,
			 std::map<std::array<int, 4>, double> &rm,
			 const double *f0,
			 double Rx, double Ry, double Rz);

double AnalyticIntegral::rep_integral(const int *pows1, const int *pows2, const int *pows3,
				      const int *pows4,
				      const std::array<double, 3> &c1,
				      const std::array<double, 3> &c2,
				      const std::array<double, 3> &c3,
				      const std::array<double, 3> &c4,
				      const GaussianOrbital *o1,
				      const GaussianOrbital *o2,
				      const GaussianOrbital *o3,
				      const GaussianOrbital *o4) const {

  double sum = 0;
  int angmom = pows1[0] + pows1[1] + pows1[2] +
    pows2[0] + pows2[1] + pows2[2] +
    pows3[0] + pows3[1] + pows3[2] +
    pows4[0] + pows4[1] + pows4[2];
  
  for(int i = 0; i < o1->getnterms(); i++) {
    for(int j = 0; j < o2->getnterms(); j++) {
      double gamma = o3->getalpha(i), delta = o4->getalpha(j);
      double eta = gamma + delta;
      double Qx = (gamma * c3[0] + delta * c4[0]) / eta,
	Qy = (gamma * c3[1] + delta * c4[1]) / eta,
	Qz = (gamma * c3[2] + delta * c4[2]) / eta;
      double Kq = M_SQRT2 * std::pow(M_PI, 1.25) / eta *
	std::exp(-gamma * delta / eta *
		 ((c3[0] - c4[0]) * (c3[0] - c4[0]) +
		  (c3[1] - c4[1]) * (c3[1] - c4[1]) +
		  (c3[2] - c4[2]) * (c3[2] - c4[2])));
      double Dq = o3->getcoef(i) * o3->getnorm(i) *
	o4->getcoef(j) * o4->getnorm(j);

      for(int k = 0; k < o1->getnterms(); k++) {
	for(int l = 0; l < o2->getnterms(); l++) {
	  double alpha = o1->getalpha(k), beta = o2->getalpha(l);
	  double zeta = alpha + beta;
	  double Px = (alpha * c1[0] + beta * c2[0]) / zeta,
	    Py = (alpha * c1[1] + beta * c2[1]) / zeta,
	    Pz = (alpha * c1[2] + beta * c2[2]) / zeta;
	  double Kp = M_SQRT2 * std::pow(M_PI, 1.25) / zeta *
	    std::exp(-alpha * beta / delta *
		     ((c1[0] - c2[0]) * (c1[0] - c2[0]) +
		      (c1[1] - c2[1]) * (c1[1] - c2[1]) +
		      (c1[2] - c2[2]) * (c1[2] - c2[2])));
	  double Dp = o1->getcoef(k) * o1->getnorm(k) *
	    o2->getcoef(l) * o2->getnorm(l);

	  double Rx = Qx - Px,
	    Ry = Qy - Py,
	    Rz = Qz - Pz;
	  double theta2 = zeta * eta / (zeta + eta);
	  double T = theta2 * (Rx * Rx + Ry * Ry + Rz * Rz);
	  double omega = Kp * Kq * Dp * Dq /
	    (std::sqrt(zeta + eta));


	  // Precompute the integrals.
	  double *f0_ints = new double[angmom + 1];
	  for(int m = 0; m <= angmom; m++) {
	    f0_ints[m] = omega * std::pow(2 * theta2, m) * boys_square(m, T);
	  }

	  std::array<int, 12> index = {
	    pows1[0], pows1[1], pows1[2],
	    pows2[0], pows2[1], pows2[2],
	    pows3[0], pows3[1], pows3[2],
	    pows4[0], pows4[1], pows4[2]};

	  std::map<std::array<int, 12>, double> abcd, apcq;
	  std::map<std::array<int, 4>, double> rm;
	  double res = compute_abcd(index, abcd, apcq, rm, f0_ints,
				    alpha, beta, gamma, delta,
				    Rx, Ry, Rz,
				    c1, c2, c3, c4);
	  sum += res;
	  delete[] f0_ints;
	}
      }
    }
  }
  return sum;
}


double AnalyticIntegral::repulsion(const GaussianOrbital *o1,
				 const GaussianOrbital *o2,
				 const GaussianOrbital *o3,
				 const GaussianOrbital *o4,
				 std::array<double, 3> c1,
				 std::array<double, 3> c2,
				 std::array<double, 3> c3,
				 std::array<double, 3> c4) const {
  double sum = 0;

  for(int i = 0; i < o1->getharms().getsize(); i++) {
    for(int j = 0; j < o2->getharms().getsize(); j++) {
      for(int k = 0; k < o3->getharms().getsize(); k++) {
	for(int l = 0; l < o4->getharms().getsize(); l++) {
	  sum += o1->getharms().getcoef(i) *
	    o2->getharms().getcoef(j) *
	    o3->getharms().getcoef(k) *
	    o4->getharms().getcoef(l) *
	    rep_integral(o1->getharms().gettermorder(i),
			 o2->getharms().gettermorder(j),
			 o3->getharms().gettermorder(k),
			 o4->getharms().gettermorder(l),
			 c1, c2, c3, c4,
			 o1, o2, o3, o4);
	}
      }
    }
  }
  return sum;
}

static double compute_abcd(const std::array<int, 12> &index,
			   std::map<std::array<int, 12>, double> &abcd,
			   std::map<std::array<int, 12>, double> &apcq,
			   std::map<std::array<int, 4>, double> &rm,
			   const double *f0,
			   double alpha, double beta,
			   double gamma, double delta,
			   double Rx, double Ry, double Rz,
			   const std::array<double, 3> &c1,
			   const std::array<double, 3> &c2,
			   const std::array<double, 3> &c3,
			   const std::array<double, 3> &c4) {
  if(abcd.count(index) != 0) {
    return abcd.at(index);
  }

  for(int i = 0; i < 12; i++) {
    if(index[i] < 0) {
      return 0;
    }
  }

  if(index[3] == 0) {
    if(index[4] == 0) {
      if(index[5] == 0) {
	if(index[9] == 0) {
	  if(index[10] == 0) {
	    if(index[11] == 0) {
	      std::array<int, 12> ind1 = {
		index[0], index[1], index[2],
		0, 0, 0,
		index[6], index[7], index[8],
		0, 0, 0};
	      double res = compute_apcq(ind1, apcq, rm, f0, alpha, beta,
					gamma, delta, Rx, Ry, Rz,
					c1, c2, c3, c4);
	      abcd[index] = res;
	      return res;
	    } else {
	      std::array<int, 12> ind1 = {
		index[0], index[1], index[2],
		index[3], index[4], index[5],
		index[6], index[7], index[8] + 1,
		index[9], index[10], index[11] - 1},
		ind2 = {
		index[0], index[1], index[2],
		index[3], index[4], index[5],
		index[6], index[7], index[8],
		index[9], index[10], index[11] - 1};
	      double res = compute_abcd(ind1, abcd, apcq, rm, f0, alpha, beta,
					gamma, delta, Rx, Ry, Rz,
					c1, c2, c3, c4) +
		(c3[2] - c4[2]) * compute_abcd(ind2, abcd, apcq, rm, f0,
					       alpha, beta,
					       gamma, delta, Rx, Ry, Rz,
					       c1, c2, c3, c4);
	      abcd[index] = res;
	      return res;
	    }
	  } else {
	    std::array<int, 12> ind1 = {
	      index[0], index[1], index[2],
	      index[3], index[4], index[5],
	      index[6], index[7] + 1, index[8],
	      index[9], index[10] - 1, index[11]},
	      ind2 = {
		index[0], index[1], index[2],
		index[3], index[4], index[5],
		index[6], index[7], index[8],
		index[9], index[10] - 1, index[11]};
	    double res = compute_abcd(ind1, abcd, apcq, rm, f0, alpha, beta,
				      gamma, delta, Rx, Ry, Rz,
				      c1, c2, c3, c4) +
	      (c3[1] - c4[1]) * compute_abcd(ind2, abcd, apcq, rm, f0,
					     alpha, beta,
					     gamma, delta, Rx, Ry, Rz,
					     c1, c2, c3, c4);
	    abcd[index] = res;
	    return res;
	  }
	} else {
	  std::array<int, 12> ind1 = {
	    index[0], index[1], index[2],
	    index[3], index[4], index[5],
	    index[6] + 1, index[7], index[8],
	    index[9] - 1, index[10], index[11]},
	    ind2 = {
	      index[0], index[1], index[2],
	      index[3], index[4], index[5],
	      index[6], index[7], index[8],
	      index[9] - 1, index[10], index[11]};
	  double res = compute_abcd(ind1, abcd, apcq, rm, f0, alpha, beta,
				    gamma, delta, Rx, Ry, Rz,
				    c1, c2, c3, c4) +
	    (c3[0] - c4[0]) * compute_abcd(ind2, abcd, apcq, rm, f0,
					   alpha, beta,
					   gamma, delta, Rx, Ry, Rz,
					   c1, c2, c3, c4);
	  abcd[index] = res;
	  return res;
	}
      } else {
	std::array<int, 12> ind1 = {
	  index[0], index[1], index[2] + 1,
	  index[3], index[4], index[5] - 1,
	  index[6], index[7], index[8],
	  index[9], index[10], index[11]},
	  ind2 = {
	    index[0], index[1], index[2],
	    index[3], index[4], index[5] - 1,
	    index[6], index[7], index[8],
	    index[9], index[10], index[11]};
	double res = compute_abcd(ind1, abcd, apcq, rm, f0, alpha, beta,
				  gamma, delta, Rx, Ry, Rz,
				  c1, c2, c3, c4) +
	  (c1[2] - c2[2]) * compute_abcd(ind2, abcd, apcq, rm, f0, alpha, beta,
					 gamma, delta, Rx, Ry, Rz,
					 c1, c2, c3, c4);
	abcd[index] = res;
	return res;
      }
    } else {
      std::array<int, 12> ind1 = {
	index[0], index[1] + 1, index[2],
	index[3], index[4] - 1, index[5],
	index[6], index[7], index[8],
	index[9], index[10], index[11]},
	ind2 = {
	  index[0], index[1], index[2],
	  index[3], index[4] - 1, index[5],
	  index[6], index[7], index[8],
	  index[9], index[10], index[11]};
      double res = compute_abcd(ind1, abcd, apcq, rm, f0, alpha, beta,
				gamma, delta, Rx, Ry, Rz,
				c1, c2, c3, c4) +
	(c1[1] - c2[1]) * compute_abcd(ind2, abcd, apcq, rm, f0, alpha, beta,
				       gamma, delta, Rx, Ry, Rz,
				       c1, c2, c3, c4);
      abcd[index] = res;
      return res;
    }
  } else {
    std::array<int, 12> ind1 = {
      index[0] + 1, index[1], index[2],
      index[3] - 1, index[4], index[5],
      index[6], index[7], index[8],
      index[9], index[10], index[11]},
      ind2 = {
	index[0], index[1], index[2],
	index[3] - 1, index[4], index[5],
	index[6], index[7], index[8],
	index[9], index[10], index[11]};
    double res = compute_abcd(ind1, abcd, apcq, rm, f0, alpha, beta,
			      gamma, delta, Rx, Ry, Rz,
			      c1, c2, c3, c4) +
      (c1[0] - c2[0]) * compute_abcd(ind2, abcd, apcq, rm, f0, alpha, beta,
				     gamma, delta, Rx, Ry, Rz,
				     c1, c2, c3, c4);
    abcd[index] = res;
    return res;
  }
}

static double compute_apcq(const std::array<int, 12> &index,
			   std::map<std::array<int, 12>, double> &apcq,
			   std::map<std::array<int, 4>, double> &rm,
			   const double *f0,
			   double alpha, double beta,
			   double gamma, double delta,
			   double Rx, double Ry, double Rz,
			   const std::array<double, 3> &c1,
			   const std::array<double, 3> &c2,
			   const std::array<double, 3> &c3,
			   const std::array<double, 3> &c4) {
  if(apcq.count(index) != 0) {
    return apcq.at(index);
  }

  for(int i = 0; i < 12; i++) {
    if(index[i] < 0) {
      return 0;
    }
  }

  if(index[0] == 0) {
    if(index[1] == 0) {
      if(index[2] == 0) {
	if(index[6] == 0) {
	  if(index[7] == 0) {
	    if(index[8] == 0) {
	      int sign = ((index[9] + index[10] + index[11]) % 2)? -1: 1;
	      std::array<int, 4> ind1 = {
		index[3] + index[9],
		index[4] + index[10],
		index[5] + index[11],
		0};
	      double res = sign * compute_rm(ind1, rm, f0, Rx, Ry, Rz);
	      apcq[index] = res;
	      return res;
	    } else {
	      std::array<int, 12> ind1 = {
		index[0], index[1], index[2],
		index[3], index[4], index[5],
		index[6], index[7], index[8] - 1,
		index[9], index[10], index[11] + 1},
		ind2 = {
		  index[0], index[1], index[2],
		  index[3], index[4], index[5],
		  index[6], index[7], index[8] - 1,
		  index[9], index[10], index[11]},
		ind3 = {
		  index[0], index[1], index[2],
		  index[3], index[4], index[5],
		  index[6], index[7], index[8] - 1,
		  index[9], index[10], index[11] - 1};
	      double res = compute_apcq(ind1, apcq, rm, f0,
					alpha, beta, gamma, delta,
					Rx, Ry, Rz, c1, c2, c3, c4) /
		(2 * (gamma + delta)) - (c3[2] - c4[2]) * (2 * delta) *
		compute_apcq(ind2, apcq, rm, f0,
			     alpha, beta, gamma, delta,
			     Rx, Ry, Rz, c1, c2, c3, c4) / (2 * (gamma + delta)) +
		index[11] * compute_apcq(ind3, apcq, rm, f0,
					 alpha, beta, gamma, delta,
					 Rx, Ry, Rz, c1, c2, c3, c4);
	      apcq[index] = res;
	      return res;
	    }
	  } else {
	    std::array<int, 12> ind1 = {
	      index[0], index[1], index[2],
	      index[3], index[4], index[5],
	      index[6], index[7] - 1, index[8],
	      index[9], index[10] + 1, index[11]},
	      ind2 = {
		index[0], index[1], index[2],
		index[3], index[4], index[5],
		index[6], index[7] - 1, index[8],
		index[9], index[10], index[11]},
	      ind3 = {
		index[0], index[1], index[2],
		index[3], index[4], index[5],
		index[6], index[7] - 1, index[8],
		index[9], index[10] - 1, index[11]};
	    double res = compute_apcq(ind1, apcq, rm, f0,
				      alpha, beta, gamma, delta,
				      Rx, Ry, Rz, c1, c2, c3, c4) /
	      (2 * (gamma + delta)) - (c3[1] - c4[1]) * (2 * delta) *
	      compute_apcq(ind2, apcq, rm, f0,
			   alpha, beta, gamma, delta,
			   Rx, Ry, Rz, c1, c2, c3, c4) / (2 * (gamma + delta)) +
	      index[10] * compute_apcq(ind3, apcq, rm, f0,
				       alpha, beta, gamma, delta,
				       Rx, Ry, Rz, c1, c2, c3, c4);
	    apcq[index] = res;
	    return res;
	  }
	} else {
	  std::array<int, 12> ind1 = {
	    index[0], index[1], index[2],
	    index[3], index[4], index[5],
	    index[6] - 1, index[7], index[8],
	    index[9] + 1, index[10], index[11]},
	    ind2 = {
	      index[0], index[1], index[2],
	      index[3], index[4], index[5],
	      index[6] - 1, index[7], index[8],
	      index[9], index[10], index[11]},
	    ind3 = {
	      index[0], index[1], index[2],
	      index[3], index[4], index[5],
	      index[6] - 1, index[7], index[8],
	      index[9] - 1, index[10], index[11]};
	  double res = compute_apcq(ind1, apcq, rm, f0,
				    alpha, beta, gamma, delta,
				    Rx, Ry, Rz, c1, c2, c3, c4) /
	    (2 * (gamma + delta)) - (c3[0] - c4[0]) * (2 * delta) *
	    compute_apcq(ind2, apcq, rm, f0,
			 alpha, beta, gamma, delta,
			 Rx, Ry, Rz, c1, c2, c3, c4) / (2 * (gamma + delta)) +
	    index[9] * compute_apcq(ind3, apcq, rm, f0,
				     alpha, beta, gamma, delta,
				     Rx, Ry, Rz, c1, c2, c3, c4);
	  apcq[index] = res;
	  return res;
	}
      } else {
	std::array<int, 12> ind1 = {
	  index[0], index[1], index[2] - 1,
	  index[3], index[4], index[5] + 1,
	  index[6], index[7], index[8],
	  index[9], index[10], index[11]},
	  ind2 = {
	    index[0], index[1], index[2] - 1,
	    index[3], index[4], index[5],
	    index[6], index[7], index[8],
	    index[9], index[10], index[11]},
	  ind3 = {
	    index[0], index[1], index[2] - 1,
	    index[3], index[4], index[5] - 1,
	    index[6], index[7], index[8],
	    index[9], index[10], index[11]};
	double res = compute_apcq(ind1, apcq, rm, f0,
				  alpha, beta, gamma, delta,
				  Rx, Ry, Rz, c1, c2, c3, c4) /
	  (2 * (alpha + beta)) - (c1[2] - c2[2]) * (2 * beta) *
	  compute_apcq(ind2, apcq, rm, f0,
		       alpha, beta, gamma, delta,
		       Rx, Ry, Rz, c1, c2, c3, c4) / (2 * (alpha + beta)) +
	  index[5] * compute_apcq(ind3, apcq, rm, f0,
				   alpha, beta, gamma, delta,
				   Rx, Ry, Rz, c1, c2, c3, c4);
	apcq[index] = res;
	return res;
      }
    } else {
      std::array<int, 12> ind1 = {
	index[0], index[1] - 1, index[2],
	index[3], index[4] + 1, index[5],
	index[6], index[7], index[8],
	index[9], index[10], index[11]},
	ind2 = {
	  index[0], index[1] - 1, index[2],
	  index[3], index[4], index[5],
	  index[6], index[7], index[8],
	  index[9], index[10], index[11]},
	ind3 = {
	  index[0], index[1] - 1, index[2],
	  index[3], index[4] - 1, index[5],
	  index[6], index[7], index[8],
	  index[9], index[10], index[11]};
      double res = compute_apcq(ind1, apcq, rm, f0,
				alpha, beta, gamma, delta,
				Rx, Ry, Rz, c1, c2, c3, c4) /
	(2 * (alpha + beta)) - (c1[1] - c2[1]) * (2 * beta) *
	compute_apcq(ind2, apcq, rm, f0,
		     alpha, beta, gamma, delta,
		     Rx, Ry, Rz, c1, c2, c3, c4) / (2 * (alpha + beta)) +
	index[4] * compute_apcq(ind3, apcq, rm, f0,
				alpha, beta, gamma, delta,
				Rx, Ry, Rz, c1, c2, c3, c4);
      apcq[index] = res;
      return res;
    }
  } else {
    std::array<int, 12> ind1 = {
      index[0] - 1, index[1], index[2],
      index[3] + 1, index[4], index[5],
      index[6], index[7], index[8],
      index[9], index[10], index[11]},
      ind2 = {
	index[0] - 1, index[1], index[2],
	index[3], index[4], index[5],
	index[6], index[7], index[8],
	index[9], index[10], index[11]},
      ind3 = {
	index[0] - 1, index[1], index[2],
	index[3] - 1, index[4], index[5],
	index[6], index[7], index[8],
	index[9], index[10], index[11]};
    double res = compute_apcq(ind1, apcq, rm, f0,
			      alpha, beta, gamma, delta,
			      Rx, Ry, Rz, c1, c2, c3, c4) /
      (2 * (alpha + beta)) - (c1[0] - c2[0]) * (2 * beta) *
      compute_apcq(ind2, apcq, rm, f0,
		   alpha, beta, gamma, delta,
		   Rx, Ry, Rz, c1, c2, c3, c4) / (2 * (alpha + beta)) +
      index[3] * compute_apcq(ind3, apcq, rm, f0,
			      alpha, beta, gamma, delta,
			      Rx, Ry, Rz, c1, c2, c3, c4);
    apcq[index] = res;
    return res;
  }
}

static double compute_rm(const std::array<int, 4> &index,
			 std::map<std::array<int, 4>, double> &rm,
			 const double *f0,
			 double Rx, double Ry, double Rz) {
  if(rm.count(index) != 0) {
    return rm.at(index);
  }

  for(int i = 0; i < 4; i++) {
    if(index[i] < 0) {
      return 0;
    }
  }

  if(index[0] == 0) {
    if(index[1] == 0) {
      if(index[2] == 0) {
	return f0[index[3]];
      } else {
	std::array<int, 4> ind1 = {
	  index[0], index[1], index[2] - 1, index[3] + 1},
	  ind2 = {
	    index[0], index[1], index[2] - 2, index[3] + 1};
	double ret = Rz * compute_rm(ind1, rm, f0, Rx, Ry, Rz) -
	  (index[2] - 1) * compute_rm(ind2, rm, f0, Rx, Ry, Rz);
	rm[index] = ret;
	return ret;
      }
    } else {
      std::array<int, 4> ind1 = {
	index[0], index[1] - 1, index[2], index[3] + 1},
	ind2 = {
	  index[0], index[1] - 2, index[2], index[3] + 1};
      double ret = Ry * compute_rm(ind1, rm, f0, Rx, Ry, Rz) -
	(index[1] - 1) * compute_rm(ind2, rm, f0, Rx, Ry, Rz);
      rm[index] = ret;
      return ret;
    }
  } else {
    std::array<int, 4> ind1 = {
      index[0] - 1, index[1], index[2], index[3] + 1},
      ind2 = {
	index[0] - 2, index[1], index[2], index[3] + 1};
    double ret = Rx * compute_rm(ind1, rm, f0, Rx, Ry, Rz) -
      (index[0] - 1) * compute_rm(ind2, rm, f0, Rx, Ry, Rz);
    rm[index] = ret;
    return ret;
  }
}
	
	
