 
#include <extramath.h>
#include "../basis_sets/basis_set.hpp"
#include "../util/atom.hpp"
#include <array>
#include <cmath>
#include <map>
#include "../util/polynomial.hpp"
#include "integrals.hpp"
#include <cstdio>

using namespace compchem;

static double compute_pCm(const std::array<int, 4> &index,
			  std::map<std::array<int, 4>, double> &ints,
			  double Rx, double Ry, double Rz,
			  double omega, double zeta) {
  if(ints.find(index) != ints.end()) {
    return ints.at(index);
  }
  if(index[0] < 0 || index[1] < 0 || index[2] < 0 || index[3] < 0) {
    return 0;
  }
  if(index[0] == 0) {
    if(index[1] == 0) {
      if(index[2] == 0) {
	if(Rx == 0 && Ry == 0 && Rz == 0) { // This is a special case.
	  double res = omega * std::pow(2 * zeta, index[3]) / (2.0 * index[3] + 1);
	  ints[index] = res;
	  return res;
	} else {
	  double res = omega * std::pow(2, index[3] - 1) /
	    std::pow(Rx * Rx + Ry * Ry + Rz * Rz, -index[3] - 0.5) /
	    std::sqrt(zeta) *
	    ingamma(index[3] + 0.5, zeta * (Rx * Rx + Ry * Ry + Rz * Rz));
	  ints[index] = res;
	  return res;
	}
      } else {
	std::array<int, 4> ind1 = {
	  index[0], index[1], index[2] - 1, index[3] + 1},
	  ind2 = {
	    index[0], index[1], index[2] - 2, index[3] + 1};
	double res = Rz * compute_pCm(ind1, ints, Rx, Ry, Rz, omega, zeta) -
	  (index[2] - 1) * compute_pCm(ind2, ints, Rx, Ry, Rz, omega, zeta);
	ints[index] = res;
	return res;
      }
    } else {
      std::array<int, 4> ind1 = {
	index[0], index[1] - 1, index[2], index[3] + 1},
	ind2 = {
	  index[0], index[1] - 2, index[2], index[3] + 1};
      double res = Rz * compute_pCm(ind1, ints, Rx, Ry, Rz, omega, zeta) -
	(index[1] - 1) * compute_pCm(ind2, ints, Rx, Ry, Rz, omega, zeta);
      ints[index] = res;
      return res;
    }
  } else {
    std::array<int, 4> ind1 = {
      index[0] - 1, index[1], index[2], index[3] + 1},
      ind2 = {
	index[0] - 2, index[1], index[2], index[3] + 1};
    double res = Rz * compute_pCm(ind1, ints, Rx, Ry, Rz, omega, zeta) -
      (index[0] - 1) * compute_pCm(ind2, ints, Rx, Ry, Rz, omega, zeta);
    ints[index] = res;
    return res;
  }
}

static double compute_apCuv(const std::array<int, 8> &index,
			  std::map<std::array<int, 8>, double> &ints,
			  const std::map<std::array<int, 5>, double> &pCuv,
			  const std::array<double, 3> &c1,
			  const std::array<double, 3> &c2) {
  if(ints.find(index) != ints.end()) {
    return ints.at(index);
  }
  if(index[0] < 0 || index[1] < 0 || index[2] < 0 || index[3] < 0 ||
     index[4] < 0 || index[5] < 0) {
    return 0;
  }
  if(index[0] == 0) {
    if(index[1] == 0) {
      if(index[2] == 0) {
	std::array<int, 5> ind2 = {
	  index[3], index[4], index[5], 0, 0};
	double res = pCuv.at(ind2);
	ints[index] = res;
	return res;
      } else {
	std::array<int, 8> ind1 = {
	  index[0], index[1], index[2] - 1,
	  index[3], index[4], index[5] - 1,
	  index[6], index[7]},
	  ind2 = {
	    index[0], index[1], index[2] - 1,
	    index[3], index[4], index[5],
	    index[6] + 1, index[7] + 1},
	  ind3 = {
	    index[0], index[1], index[2] - 1,
	    index[3], index[4], index[5] + 1,
	    index[6], index[7] + 1};
	double res = index[5] * compute_apCuv(ind1, ints, pCuv, c1, c2) -
	  (c1[2] - c2[2]) * compute_apCuv(ind2, ints, pCuv, c1, c2) +
	  compute_apCuv(ind3, ints, pCuv, c1, c2);
	ints[index] = res;
	return res;
      }
    } else {
      std::array<int, 8> ind1 = {
	index[0], index[1] - 1, index[2],
	index[3], index[4] - 1, index[5],
	index[6], index[7]},
	ind2 = {
	  index[0], index[1] - 1, index[2],
	  index[3], index[4], index[5],
	  index[6] + 1, index[7] + 1},
	ind3 = {
	  index[0], index[1] - 1, index[2],
	  index[3], index[4] + 1, index[5],
	  index[6], index[7] + 1};
      double res = index[5] * compute_apCuv(ind1, ints, pCuv, c1, c2) -
	(c1[1] - c2[1]) * compute_apCuv(ind2, ints, pCuv, c1, c2) +
	compute_apCuv(ind3, ints, pCuv, c1, c2);
      ints[index] = res;
      return res;
    }
  } else {
    std::array<int, 8> ind1 = {
      index[0] - 1, index[1], index[2],
      index[3] - 1, index[4], index[5],
      index[6], index[7]},
      ind2 = {
	index[0] - 1, index[1], index[2],
	index[3], index[4], index[5],
	index[6] + 1, index[7] + 1},
      ind3 = {
	index[0] - 1, index[1], index[2],
	index[3] + 1, index[4], index[5],
	index[6], index[7] + 1};
    double res = index[5] * compute_apCuv(ind1, ints, pCuv, c1, c2) -
      (c1[0] - c2[0]) * compute_apCuv(ind2, ints, pCuv, c1, c2) +
      compute_apCuv(ind3, ints, pCuv, c1, c2);
    ints[index] = res;
    return res;
  }
}
	
static double compute_abC(const std::array<int, 6> &index,
			  std::map<std::array<int, 6>, double> &ints,
			  const std::map<std::array<int, 8>, double> &apCuv,
			  const std::array<double, 3> &c1,
			  const std::array<double, 3> &c2) {
  if(ints.find(index) != ints.end()) {
    return ints.at(index);
  }
  if(index[0] < 0 || index[1] < 0 || index[2] < 0 ||
     index[3] < 0 || index[4] < 0 || index[5] < 0) {
    return 0;
  }
  if(index[3] == 0) {
    if(index[4] == 0) {
      if(index[5] == 0) {
	std::array<int, 8> ind1 = {
	  index[0], index[1], index[2],
	  index[3], index[4], index[5],
	  0, 0};
	double res = apCuv.at(ind1);
	ints[index] = res;
	return res;
      } else {
	std::array<int, 6> ind1 = {
	  index[0], index[1], index[2] + 1,
	  index[3], index[4], index[5] - 1},
	  ind2 = {
	    index[0], index[1], index[2],
	    index[3], index[4], index[5] - 1};
	double res = compute_abC(ind1, ints, apCuv, c1, c2) +
	  (c1[2] - c2[2]) * compute_abC(ind2, ints, apCuv, c1, c2);
	ints[index] = res;
	return res;
      }
    } else {
      std::array<int, 6> ind1 = {
	index[0], index[1] + 1, index[2],
	index[3], index[4] - 1, index[5]},
	ind2 = {
	  index[0], index[1], index[2],
	  index[3], index[4] - 1, index[5]};
	double res = compute_abC(ind1, ints, apCuv, c1, c2) +
	  (c1[1] - c2[1]) * compute_abC(ind2, ints, apCuv, c1, c2);
	ints[index] = res;
	return res;
    }
  } else {
    std::array<int, 6> ind1 = {
      index[0] + 1, index[1], index[2],
      index[3] - 1, index[4], index[5]},
      ind2 = {
	index[0], index[1], index[2],
	index[3] - 1, index[4], index[5]};
    double res = compute_abC(ind1, ints, apCuv, c1, c2) +
      (c1[0] - c2[0]) * compute_abC(ind2, ints, apCuv, c1, c2);
    ints[index] = res;
    return res;
  }
}
			 

static double coul_integral(const int *pows1,
			    const int *pows2,
			    const std::array<double, 3> &c1,
			    const std::array<double, 3> &c2,
			    const GaussianOrbital *o1,
			    const GaussianOrbital *o2,
			    const Atom &atom) {
  
  std::map<std::array<int, 5>, double> pCuv;
  for(int i = 0; i < o1->getnterms(); i++) {
    for(int j = 0; j < o2->getnterms(); j++) {
      // Set up constants.
      double zeta = o1->getalpha(i) + o2->getalpha(j),
	K = 2 * M_PI * std::exp(-o1->getalpha(i) * o2->getalpha(j) / zeta *
				((c1[0] - c2[0]) * (c1[0] - c2[0]) +
				 (c1[1] - c2[1]) * (c1[1] - c2[1]) +
				 (c1[2] - c2[2]) * (c1[2] - c2[2]))) / zeta,
	Px = (o1->getalpha(i) * c1[0] + o2->getalpha(j) * c2[0]) / zeta,
	Py = (o1->getalpha(i) * c1[1] + o2->getalpha(j) * c2[1]) / zeta,
	Pz = (o1->getalpha(i) * c1[2] + o2->getalpha(j) * c2[2]) / zeta,
	Rx = atom.getx() - Px,
	Ry = atom.gety() - Py,
	Rz = atom.getz() - Pz,
	T = zeta * (Rx * Rx + Ry * Ry + Rz * Rz),
	omega = K / std::pow(2 * zeta, pows1[0] + pows1[1] + pows1[2] +
			     pows2[0] + pows2[1] + pows2[2]);
	
      // Compute [p|C](m)
      std::map<std::array<int, 4>, double> pCm;
      for(int x = 0; x <= pows1[0] + pows2[0]; x++) {
	for(int y = 0; y <= pows1[1] + pows2[1]; y++) {
	  for(int z = 0; z <= pows1[2] + pows2[2]; z++) {
	    std::array<int, 4> index = {
	      x, y, z, 0};
	    double res = compute_pCm(index, *&pCm, Rx, Ry, Rz, omega, zeta);
	    // Compute (p|C)_u,v
	    for(int u = 0; u <= pows1[0] + pows1[1] + pows1[2] +
		  pows2[0] + pows2[1] + pows2[2]; u++) {
	      for(int v = 0; v <= pows1[0] + pows1[1] + pows1[2] +
		  pows2[0] + pows2[1] + pows2[2]; v++) {
		std::array<int, 5> index2 = {
		  x, y, z, u, v};
		if(pCuv.find(index2) == pCuv.end()) {
		  pCuv[index2] = 0;
		}
		pCuv.at(index2) += o1->getcoef(i) * o2->getcoef(j) *
		  o1->getnorm(i) * o2->getnorm(j) *
		  std::pow(2 * o2->getalpha(j), u) /
		  std::pow(2 * zeta, v) * res;
	      }
	    }
	  }
	}
      }
    }
  }
  // Compute (e0|C)
  std::map<std::array<int, 8>, double> apCuv;
  for(int x = 0; x <= pows1[0] + pows2[0]; x++) {
    for(int y = 0; y <= pows1[1] + pows2[1]; y++) {
      for(int z = 0; z <= pows1[2] + pows2[2]; z++) {
	std::array<int, 8> index = {
	  x, y, z,
	  0, 0, 0,
	  0, 0};
	double res = compute_apCuv(index, *&apCuv, pCuv, c1, c2);
      }
    }
  }
  
  // Compute (ab|C)
  std::map<std::array<int, 6>, double> abC;
  std::array<int, 6> index = {
    pows1[0], pows1[1], pows1[2],
    pows2[0], pows2[1], pows2[2]};
  return compute_abC(index, *&abC, apCuv, c1, c2);
}
  
			    

double AnalyticIntegral::coulomb(const GaussianOrbital *o1,
				 const GaussianOrbital *o2,
				 std::array<double, 3> c1,
				 std::array<double, 3> c2,
				 const Atom &atom) const {
  double sum = 0;
  for(int i = 0; i < o1->getharms().getsize(); i++) {
    for(int j = 0; j < o2->getharms().getsize(); j++) {
      sum += o1->getharms().getcoef(i) * o2->getharms().getcoef(j) *
	coul_integral(o1->getharms().gettermorder(i),
		      o2->getharms().gettermorder(j),
		      c1, c2, o1, o2, atom);
    }
  }
  return -sum * atom.getZ();
}
