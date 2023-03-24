
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

static double compute_ap(const std::array<int, 6> &index,
			 std::map<std::array<int, 6>, double> &ints,
			 double zeta,
			 double Ax, double Ay, double Az,
			 double Px, double Py, double Pz,
			 double K) {
  if(ints.find(index) != ints.end()) {
    return ints.at(index);
  }
  if(index[0] < 0 || index[1] < 0 || index[2] < 0 ||
     index[3] < 0 || index[4] < 0 || index[5] < 0) {
    return 0;
  }
  if(index[0] == 0) {
    if(index[1] == 0) {
      if(index[2] == 0) {
	if(index[3] % 2 == 1 || index[4] % 2 == 1 || index[5] % 2 == 1) {
	  ints[index] = 0;
	  return 0;
	}
	double res = K * std::tgamma(index[3] / 2 + 0.5) *
	  std::tgamma(index[4] / 2 + 0.5) *
	  std::tgamma(index[5] / 2 + 0.5) /
	  std::pow(zeta, (index[3] + index[4] + index[5]) / 2 + 1.5);
	ints[index] = res;
	return res;
      } else {
	std::array<int, 6> ind1 = {
	  index[0], index[1], index[2] - 1,
	  index[3], index[4], index[5] + 1},
	  ind2 = {
	    index[0], index[1], index[2] - 1,
	    index[3], index[4], index[5]};
	double res = compute_ap(ind1, ints, zeta, Ax, Ay, Az, Px, Py, Pz, K) +
	  (Pz - Az) * compute_ap(ind2, ints, zeta, Ax, Ay, Az, Px, Py, Pz, K);
	ints[index] = res;
	return res;
      }
    } else {
      std::array<int, 6> ind1 = {
	index[0], index[1] - 1, index[2],
	index[3], index[4] + 1, index[5]},
	ind2 = {
	  index[0], index[1] - 1, index[2],
	  index[3], index[4], index[5]};
      double res = compute_ap(ind1, ints, zeta, Ax, Ay, Az, Px, Py, Pz, K) +
	(Py - Ay) * compute_ap(ind2, ints, zeta, Ax, Ay, Az, Px, Py, Pz, K);
      ints[index] = res;
      return res;
    }
  } else {
    std::array<int, 6> ind1 = {
      index[0] - 1, index[1], index[2],
      index[3] + 1, index[4], index[5]},
      ind2 = {
	index[0] - 1, index[1], index[2],
	index[3], index[4], index[5]};
    double res = compute_ap(ind1, ints, zeta, Ax, Ay, Az, Px, Py, Pz, K) +
      (Px - Ax) * compute_ap(ind2, ints, zeta, Ax, Ay, Az, Px, Py, Pz, K);
    ints[index] = res;
    return res;
  }
}

static double compute_ab(const std::array<int, 6> &index,
			 std::map<std::array<int, 6>, double> &ints,
			 const std::map<std::array<int, 3>, double> &e0con,
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
	std::array<int, 3> ind = {
	  index[0], index[1], index[2]};
	double res = e0con.at(ind);
	ints[index] = res;
	return res;
      } else {
	std::array<int, 6> ind1 = {
	  index[0], index[1], index[2] + 1,
	  index[3], index[4], index[5] - 1},
	  ind2 = {
	    index[0], index[1], index[2],
	    index[3], index[4], index[5] - 1};
	double res = compute_ab(ind1, ints, e0con, c1, c2) +
	  (c1[2] - c2[2]) * compute_ab(ind2, ints, e0con, c1, c2);
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
      double res = compute_ab(ind1, ints, e0con, c1, c2) +
	(c1[1] - c2[1]) * compute_ab(ind2, ints, e0con, c1, c2);
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
    double res = compute_ab(ind1, ints, e0con, c1, c2) +
      (c1[0] - c2[0]) * compute_ab(ind2, ints, e0con, c1, c2);
    ints[index] = res;
    return res;
  }
}

static double compute_overlap(const int *pow1, const int *pow2,
			      const GaussianOrbital *o1,
			      const GaussianOrbital *o2,
			      const std::array<double, 3> &center1,
			      const std::array<double, 3> &center2) {
  std::map<std::array<int, 3>, double> e0con;
  for(int i = 0; i < o1->getnterms(); i++) {
    for(int j = 0; j < o2->getnterms(); j++) {
      double zeta = o1->getalpha(i) + o2->getalpha(j),
	Px = (o1->getalpha(i) * center1[0] + o2->getalpha(j) * center2[0]) / zeta,
	Py = (o1->getalpha(i) * center1[1] + o2->getalpha(j) * center2[1]) / zeta,
	Pz = (o1->getalpha(i) * center1[2] + o2->getalpha(j) * center2[2]) / zeta,
	K = o1->getcoef(i) * o1->getnorm(i) * o2->getcoef(j) * o2->getnorm(j) *
	std::exp(-o1->getalpha(i) * o2->getalpha(j) / zeta *
		 ((center1[0] - center2[0]) * (center1[0] - center2[0]) +
		  (center1[1] - center2[1]) * (center1[1] - center2[1]) +
		  (center1[2] - center2[2]) * (center1[2] - center2[2])));
      // Compute the uncontracted integrals [e|0].
      std::map<std::array<int, 6>, double> ap;
      for(int k = 0; k <= pow1[0] + pow2[0]; k++) {
	for(int l = 0; l <= pow1[1] + pow2[1]; l++) {
	  for(int m = 0; m <= pow1[2] + pow2[2]; m++) {
	    std::array<int, 6> index = {k, l, m, 0, 0, 0};
	    double res = compute_ap(*&index, *&ap, zeta,
				    center1[0], center1[1], center1[2],
				    Px, Py, Pz, K);
	    // Compute the contracted integrals (e|0).
	    std::array<int, 3> index2 = {k, l, m};
	    if(e0con.find(index2) == e0con.end()) {
	      e0con[index2] = 0;
	    }
	    e0con.at(index2) += res;
	  }
	}
      }
    }
  }
  std::map<std::array<int, 6>, double> ab;

  std::array<int, 6> index = {
    pow1[0], pow1[1], pow1[2],
    pow2[0], pow2[1], pow2[2]};
  return compute_ab(*&index, *&ab, *&e0con, center1, center2);
}
	

double AnalyticIntegral::overlap(const GaussianOrbital *o1,
				 const GaussianOrbital *o2,
				 std::array<double, 3> center1,
				 std::array<double, 3> center2) const {
  double sum = 0;
  for(int i = 0; i < o1->getharms().getsize(); i++) {
    for(int j = 0; j < o2->getharms().getsize(); j++) {
      sum += o1->getharms().getcoef(i) * o2->getharms().getcoef(j) *
	compute_overlap(o1->getharms().gettermorder(i),
			o2->getharms().gettermorder(j),
			o1, o2, center1, center2);
    }
  }
  return sum;
  
}
