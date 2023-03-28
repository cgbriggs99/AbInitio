 
#include <extramath.h>
#include "../basis_sets/basis_set.hpp"
#include "../util/atom.hpp"
#include <array>
#include <cmath>
#include <vector>
#include <map>
#include "../util/polynomial.hpp"
#include "integrals.hpp"
#include <cstdio>

using namespace compchem;    

double AnalyticIntegral::compute_rm(const std::array<int, 4> &index,
			 std::map<std::array<int, 4>, double> &ints,
			 double Rx, double Ry, double Rz, double omega,
			 double theta2, double T) const {
  if(ints.count(index) == 1) {
    return ints.at(index);
  }
  if(index[0] < 0 || index[1] < 0 || index[2] < 0 || index[3] < 0) {
    return 0;
  }

  if(index[0] == 0) {
    if(index[1] == 0) {
      if(index[2] == 0) {
	if(T == 0) {
	  double res = omega * std::pow(2 * theta2, index[3]) /
	    (2 * index[3] + 1);
	  ints[index] = res;
	  return res;
	} else {
	  double res = omega * std::pow(2 * theta2, index[3]) *
	    boys_square(index[3], T);
	  ints[index] = res;
	  return res;
	}
      } else {
	std::array<int, 4> ind1 = {
	  index[0], index[1], index[2] - 1, index[3] + 1},
	  ind2 = {
	    index[0], index[1], index[2] - 2, index[3] + 1};
	double res = Rz * compute_rm(ind1, ints, Rx, Ry, Rz, omega, theta2, T) -
	  (index[2] - 1) * compute_rm(ind2, ints, Rx, Ry, Rz, omega, theta2, T);
	ints[index] = res;
	return res;
      }
    } else {
      std::array<int, 4> ind1 = {
	index[0], index[1] - 1, index[2], index[3] + 1},
	ind2 = {
	  index[0], index[1] - 2, index[2], index[3] + 1};
      double res = Ry * compute_rm(ind1, ints, Rx, Ry, Rz, omega, theta2, T) -
	(index[1] - 1) * compute_rm(ind2, ints, Rx, Ry, Rz, omega, theta2, T);
      ints[index] = res;
      return res;
    }
  } else {
    std::array<int, 4> ind1 = {
      index[0] - 1, index[1], index[2], index[3] + 1},
      ind2 = {
	index[0] - 2, index[1], index[2], index[3] + 1};
    double res = Rx * compute_rm(ind1, ints, Rx, Ry, Rz, omega, theta2, T) -
      (index[0] - 1) * compute_rm(ind2, ints, Rx, Ry, Rz, omega, theta2, T);
    ints[index] = res;
    return res;
  }
}

double AnalyticIntegral::compute_apquv(const std::array<int, 11> &index,
			    std::map<std::array<int, 11>, double> &ints,
			    const std::map<std::array<int, 8>, double> &pquv,
			    const std::array<double, 3> &c1,
			    const std::array<double, 3> &c2) const {
  if(ints.count(index) != 0) {
    return ints.at(index);
  }
  for(int i = 0; i < 11; i++) {
    if(index[i] < 0) {
      return 0;
    }
  }

  if(index[0] == 0) {
    if(index[1] == 0) {
      if(index[2] == 0) {
	std::array<int, 8> ind1 = {
	  index[3], index[4], index[5],
	  index[6], index[7], index[8],
	  index[9], index[10]};
	ints[index] = pquv.at(ind1);
	return pquv.at(ind1);
      } else {
	std::array<int, 11> ind1 = {
	  index[0], index[1], index[2] - 1,
	  index[3], index[4], index[5] - 1,
	  index[6], index[7], index[8],
	  index[9], index[10]},
	  ind2 = {
	    index[0], index[1], index[2] - 1,
	    index[3], index[4], index[5],
	    index[6], index[7], index[8],
	    index[9] + 1, index[10] + 1},
	  ind3 = {
	    index[0], index[1], index[2] - 1,
	    index[3], index[4], index[5] + 1,
	    index[6], index[7], index[8],
	    index[9], index[10] + 1};
	double res = index[2] * compute_apquv(ind1, ints, pquv, c1, c2) -
	  (c1[2] - c2[2]) * compute_apquv(ind2, ints, pquv, c1, c2) +
	  compute_apquv(ind3, ints, pquv, c1, c2);
	ints[index] = res;
	return res;
      }
    } else {
      std::array<int, 11> ind1 = {
	index[0], index[1] - 1, index[2],
	index[3], index[4] - 1, index[5],
	index[6], index[7], index[8],
	index[9], index[10]},
	ind2 = {
	  index[0], index[1] - 1, index[2],
	  index[3], index[4], index[5],
	  index[6], index[7], index[8],
	  index[9] + 1, index[10] + 1},
	ind3 = {
	  index[0], index[1] - 1, index[2],
	  index[3], index[4] + 1, index[5],
	  index[6], index[7], index[8],
	  index[9], index[10] + 1};
      double res = index[1] * compute_apquv(ind1, ints, pquv, c1, c2) -
	(c1[1] - c2[1]) * compute_apquv(ind2, ints, pquv, c1, c2) +
	compute_apquv(ind3, ints, pquv, c1, c2);
      ints[index] = res;
      return res;
    }
  } else {
    std::array<int, 11> ind1 = {
      index[0] - 1, index[1], index[2],
      index[3] - 1, index[4], index[5],
      index[6], index[7], index[8],
      index[9], index[10]},
      ind2 = {
	index[0] - 1, index[1], index[2],
	index[3], index[4], index[5],
	index[6], index[7], index[8],
	index[9] + 1, index[10] + 1},
      ind3 = {
	index[0] - 1, index[1], index[2],
	index[3] + 1, index[4], index[5],
	index[6], index[7], index[8],
	index[9], index[10] + 1};
    double res = index[0] * compute_apquv(ind1, ints, pquv, c1, c2) -
      (c1[0] - c2[0]) * compute_apquv(ind2, ints, pquv, c1, c2) +
      compute_apquv(ind3, ints, pquv, c1, c2);
    ints[index] = res;
    return res;
  }
}

double AnalyticIntegral::compute_e0cquv(const std::array<int, 11> &index,
			     std::map<std::array<int, 11>, double> &ints,
			     const std::map<std::array<int, 8>, double> &e0quv,
			     const std::array<double, 3> &c3,
			     const std::array<double, 3> &c4) const {
  if(ints.count(index) != 0) {
    return ints[index];
  }

  for(int i = 0; i < 11; i++) {
    if(index[i] < 0) {
      return 0;
    }
  }

  if(index[3] == 0) {
    if(index[4] == 0) {
      if(index[5] == 0) {
	std::array<int, 8> ind2 = {
	  index[0], index[1], index[2],
	  index[6], index[7], index[8],
	  index[9], index[10]};
	return e0quv.at(ind2);
      } else {
	std::array<int, 11> ind1 = {
	  index[0], index[1], index[2] - 1,
	  index[3], index[4], index[5] - 1,
	  index[6], index[7], index[8],
	  index[9], index[10]},
	  ind2 = {
	  index[0], index[1], index[2] - 1,
	  index[3], index[4], index[5],
	  index[6], index[7], index[8],
	  index[9] + 1, index[10] + 1},
	  ind3 = {
	  index[0], index[1], index[2] - 1,
	  index[3], index[4], index[5] + 1,
	  index[6], index[7], index[8],
	  index[9], index[10] + 1};
	double res = index[2] * compute_e0cquv(ind1, ints, e0quv, c3, c4) -
	  (c3[2] - c4[2]) * compute_e0cquv(ind2, ints, e0quv, c3, c4) +
	  compute_e0cquv(ind3, ints, e0quv, c3, c4);
	ints[index] = res;
	return res;
      }
    } else {
      std::array<int, 11> ind1 = {
	index[0], index[1] - 1, index[2],
	index[3], index[4] - 1, index[5],
	index[6], index[7], index[8],
	index[9], index[10]},
	ind2 = {
	  index[0], index[1] - 1, index[2],
	  index[3], index[4], index[5],
	  index[6], index[7], index[8],
	  index[9] + 1, index[10] + 1},
	ind3 = {
	  index[0], index[1] - 1, index[2],
	  index[3], index[4] + 1, index[5],
	  index[6], index[7], index[8],
	  index[9], index[10] + 1};
      double res = index[1] * compute_e0cquv(ind1, ints, e0quv, c3, c4) -
	(c3[1] - c4[1]) * compute_e0cquv(ind2, ints, e0quv, c3, c4) +
	compute_e0cquv(ind3, ints, e0quv, c3, c4);
      ints[index] = res;
      return res;
    }
  } else {
    std::array<int, 11> ind1 = {
      index[0] - 1, index[1], index[2],
      index[3] - 1, index[4], index[5],
      index[6], index[7], index[8],
      index[9], index[10]},
      ind2 = {
	index[0] - 1, index[1], index[2],
	index[3], index[4], index[5],
	index[6], index[7], index[8],
	index[9] + 1, index[10] + 1},
      ind3 = {
	index[0] - 1, index[1], index[2],
	index[3] + 1, index[4], index[5],
	index[6], index[7], index[8],
	index[9], index[10] + 1};
    double res = index[0] * compute_e0cquv(ind1, ints, e0quv, c3, c4) -
      (c3[0] - c4[0]) * compute_e0cquv(ind2, ints, e0quv, c3, c4) +
      compute_e0cquv(ind3, ints, e0quv, c3, c4);
    ints[index] = res;
    return res;
  }
}

double AnalyticIntegral::compute_abcd(const std::array<int, 12> &index,
			   std::map<std::array<int, 12>, double> &ints,
			   const std::map<std::array<int, 6>, double> &e0f0,
			   const std::array<double, 3> &c1,
			   const std::array<double, 3> &c2,
			   const std::array<double, 3> &c3,
			   const std::array<double, 3> &c4) const {
  if(ints.count(index) != 0) {
    return ints.at(index);
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
	      std::array<int, 6> ind = {
		index[0], index[1], index[2],
		index[6], index[7], index[8]};
	      double res = e0f0.at(ind);
	      ints[index] = res;
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
	      double res = compute_abcd(ind1, ints, e0f0, c1, c2, c3, c4) +
		(c3[2] - c4[2]) *
		compute_abcd(ind2, ints, e0f0, c1, c2, c3, c4);
	      ints[index] = res;
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
	      double res = compute_abcd(ind1, ints, e0f0, c1, c2, c3, c4) +
		(c3[1] - c4[1]) *
		compute_abcd(ind2, ints, e0f0, c1, c2, c3, c4);
	      ints[index] = res;
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
	  double res = compute_abcd(ind1, ints, e0f0, c1, c2, c3, c4) +
	    (c3[0] - c4[0]) *
	    compute_abcd(ind2, ints, e0f0, c1, c2, c3, c4);
	  ints[index] = res;
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
	double res = compute_abcd(ind1, ints, e0f0, c1, c2, c3, c4) +
	  (c1[2] - c2[2]) *
	  compute_abcd(ind2, ints, e0f0, c1, c2, c3, c4);
	ints[index] = res;
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
      double res = compute_abcd(ind1, ints, e0f0, c1, c2, c3, c4) +
	(c1[1] - c2[1]) *
	compute_abcd(ind2, ints, e0f0, c1, c2, c3, c4);
      ints[index] = res;
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
    double res = compute_abcd(ind1, ints, e0f0, c1, c2, c3, c4) +
      (c1[0] - c2[0]) *
      compute_abcd(ind2, ints, e0f0, c1, c2, c3, c4);
    ints[index] = res;
    return res;
  }
}

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
  std::map<std::array<int, 8>, double> e0quv;
  for(int i = 0; i < o3->getnterms(); i++) {
    for(int j = 0; j < o4->getnterms(); j++) {
      // Define constants.
      double eta = o3->getalpha(i) + o4->getalpha(j);
      double Qx = (o3->getalpha(i) * c3[0] + o4->getalpha(j) * c4[0]) / eta,
	Qy = (o3->getalpha(i) * c3[1] + o4->getalpha(j) * c4[1]) / eta,
	Qz = (o3->getalpha(i) * c3[2] + o4->getalpha(j) * c4[2]) / eta;
      double Kq = M_SQRT2 * std::pow(M_PI, 1.25) / eta *
	std::exp(-o3->getalpha(i) * o4->getalpha(j) / eta *
		 ((c3[0] - c4[0]) * (c3[0] - c4[0]) +
		  (c3[1] - c4[1]) * (c3[1] - c4[1]) +
		  (c3[2] - c4[2]) * (c3[2] - c4[2])));
      double Dq = o3->getcoef(i) * o4->getcoef(j) *
	o3->getnorm(i) * o4->getnorm(j);

      std::map<std::array<int, 8>, double> pquv;

      for(int k = 0; k < o1->getnterms(); k++) {
	for(int l = 0; l < o2->getnterms(); l++) {
	  // More constants
	  double zeta = o1->getalpha(k) + o2->getalpha(l);
	  double Px =
	    (o1->getalpha(k) * c1[0] + o2->getalpha(l) * c2[0]) / zeta,
	    Py = (o1->getalpha(k) * c1[1] + o2->getalpha(l) * c2[1]) / zeta,
	    Pz = (o1->getalpha(k) * c1[2] + o2->getalpha(l) * c2[2]) / zeta;
	  double Kp = M_SQRT2 * std::pow(M_PI, 1.25) / zeta *
	    std::exp(-o1->getalpha(k) * o2->getalpha(l) / zeta *
		     ((c1[0] - c2[0]) * (c1[0] - c2[0]) +
		      (c1[1] - c2[1]) * (c1[1] - c2[1]) +
		      (c1[2] - c2[2]) * (c1[2] - c2[2])));
	  double Dp = o1->getcoef(k) * o2->getcoef(l) *
	    o1->getnorm(k) * o2->getnorm(l);
	  double Rx = Qx - Px,
	    Ry = Qy - Py,
	    Rz = Qz - Pz;
	  double theta2 = zeta * eta / (zeta + eta);
	  double T = theta2 * (Rx * Rx + Ry * Ry + Rz * Rz);
	  double omega = Kp * Dp * Kq * Dq /
	    (std::pow(2 * zeta, pows1[0] + pows1[1] + pows1[2] +
		      pows2[0] + pows2[1] + pows2[2]) *
	     std::pow(2 * eta, pows3[0] + pows3[1] + pows3[2] +
		      pows4[0] + pows4[1] + pows4[2]) *
	     std::sqrt(zeta + eta));
	    
	  // Generate [r](m).
	  std::map<std::array<int, 4>, double> rm;
	  for(int px = 0; px <= pows1[0] + pows2[0] + pows3[0] + pows4[0]; px++) {
	    for(int py = 0; py <= pows1[1] + pows2[1] + pows3[1] + pows4[1]; py++) {
	      for(int pz = 0; pz <= pows1[2] + pows2[2] + pows3[2] + pows4[2]; pz++) {
		for(int m = 0; m <= o1->getl() + o2->getl() + o3->getl() + o4->getl(); m++) {
		  std::array<int, 4> index = {
		    px, py, pz, m};
		  double res = compute_rm(index, *&rm,
					  Rx, Ry, Rz, omega, theta2, T);
		}
	      }
	    }
	  }
	  
	  // Generate (p|q]uv. May be optimized for storage somehow.
	  for(int px = 0; px <= pows1[0] + pows2[0]; px++) {
	    for(int py = 0; py <= pows1[1] + pows2[1]; py++) {
	      for(int pz = 0; pz <= pows1[2] + pows2[2]; pz++) {
		for(int qx = 0; qx <= pows3[0] + pows4[0]; qx++) {
		  for(int qy = 0; qy <= pows3[1] + pows4[1]; qy++) {
		    for(int qz = 0; qz <= pows3[2] + pows4[2]; qz++) {
		      for(int u = 0; u <= pows1[0] + pows1[1] + pows1[2] +
			    pows2[0] + pows2[1] + pows2[2]; u++) {
			for(int v = 0; v <= pows1[0] + pows1[1] + pows1[2] +
			      pows2[0] + pows2[1] + pows2[2]; v++) {
			  std::array<int, 8> index = {
			    px, py, pz, qx, qy, qz, u, v};
			  std::array<int, 4> index2 = {
			    px + qx, py + qy, pz + qz, 0};
			  if(pquv.count(index) == 0) {
			    pquv[index] = 0;
			  }
			  if((qx + qy + qz) % 2 == 1) {
			    pquv.at(index) -=
			      std::pow(2 * o2->getalpha(l), u) /
			      std::pow(2 * zeta, v) * rm.at(index2);
			  } else {
			    pquv.at(index) +=
			      std::pow(2 * o2->getalpha(l), u) /
			      std::pow(2 * zeta, v) * rm.at(index2);
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }

      // Generate (e0|q].
      std::map<std::array<int, 6>, double> e0q;
      std::map<std::array<int, 11>, double> apquv;

      for(int ex = 0; ex <= pows1[0] + pows2[0]; ex++) {
	for(int ey = 0; ey <= pows1[1] + pows2[1]; ey++) {
	  for(int ez = 0; ez <= pows1[2] + pows2[2]; ez++) {
	    for(int qx = 0; qx <= pows3[0] + pows4[0]; qx++) {
	      for(int qy = 0; qy <= pows3[1] + pows4[1]; qy++) {
		for(int qz = 0; qz <= pows3[2] + pows4[2]; qz++) {
		  std::array<int, 11> ind1 = {
		    ex, ey, ez, 0, 0, 0, qx, qy, qz, 0, 0};
		  double res = compute_apquv(*&ind1, *&apquv, *&pquv, c1, c2);
		  std::array<int, 6> ind2 = {
		    ex, ey, ez, qx, qy, qz};
		  e0q[ind2] = res;

		  // Contract to (e0|q)_uv.
		  for(int u = 0; u <= pows3[0] + pows3[1] + pows3[2] +
			pows4[0] + pows4[1] + pows4[2]; u++) {
		    for(int v = 0; v <= pows3[0] + pows3[1] + pows3[2] +
			  pows4[0] + pows4[1] + pows4[2]; v++) {
		      std::array<int, 8> index = {ex, ey, ez, qx, qy, qz, u, v};
		      
		      if(e0quv.count(index) == 0) {
			e0quv[index] = 0;
		      }
		      e0quv[index] = e0quv.at(index) +
			std::pow(2 * o4->getalpha(j), u) /
			std::pow(2 * eta, v) * res;
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }

  // Generate (e0|f0).
  std::map<std::array<int, 6>, double> e0f0;
  std::map<std::array<int, 11>, double> e0cquv;
  
  for(int ex = 0; ex <= pows1[0] + pows2[0]; ex++) {
    for(int ey = 0; ey <= pows1[1] + pows2[1]; ey++) {
      for(int ez = 0; ez <= pows1[2] + pows2[2]; ez++) {
	for(int fx = 0; fx <= pows3[0] + pows4[0]; fx++) {
	  for(int fy = 0; fy <= pows3[1] + pows4[1]; fy++) {
	    for(int fz = 0; fz <= pows3[2] + pows4[2]; fz++) {
	      std::array<int, 11> index = {
		ex, ey, ez, fx, fy, fz, 0, 0, 0, 0, 0};
	      double res =
		compute_e0cquv(*&index, *&e0cquv, *&e0quv, c3, c4);
	      std::array<int, 6> ind2 = {
		ex, ey, ez, fx, fy, fz};
	      e0f0[ind2] = res;
	    }
	  }
	}
      }
    }
  }

  // Generate (ab|cd).
  std::map<std::array<int, 12>, double> abcd;
  std::array<int, 12> index = {
    pows1[0], pows1[1], pows1[2],
    pows2[0], pows2[1], pows2[2],
    pows3[0], pows3[1], pows3[2],
    pows4[0], pows4[1], pows4[2]};
  return compute_abcd(*&index, *&abcd, *&e0f0, c1, c2, c3, c4);
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
