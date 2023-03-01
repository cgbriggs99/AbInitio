 
#include <extramath.h>
#include "../basis_sets/basis_sets.hpp"
#include "../util/atom.hpp"
#include <array>
#include <vector>
#include <map>
#include <initializer_list>
#include <cmath>
#include "../util/polynomial.hpp"
#include "integrals.hpp"


using namespace compchem;

static double compute_pq(std::array<int, 4> pqm,
			 std::map<std::array<int, 4>, double> &ints,
			 double omega, double theta2,
			 double T, double Rx, double Ry, double Rz) {
  if(ints.find(pqm) != ints.end()) {
    return ints.at(pqm);
  }
  if(pqm[0] < 0 || pqm[1] < 0 || pqm[2] < 0) {
    return 0;
  }
  if(pqm[0] == 0 && pqm[1] == 0 && pqm[2] == 0) {
    // If non-zero, incomplete gamma function.
    if(T != 0) {
      return omega * std::pow(2 * theta2, pqm[3]) *
	std::pow(T, -2 * pqm[3] - 1) * ingamma(2 * pqm[3] + 1, 1 / T);
    } else {  // If zero, simple monomial.
      return omega * std::pow(2 * theta2, pqm[3]) / (2 * pqm[3] + 1);
    }
  } else {
    double res = Rx *
      compute_pq(std::array<int, 4>({pqm[0] - 1, pqm[1], pqm[2], pqm[3] + 1}),
		 ints, omega, theta2, T, Rx, Ry, Rz) +
      Ry * compute_pq(std::array<int, 4>
		      ({pqm[0], pqm[1] - 1, pqm[2], pqm[3] + 1}),
		      ints, omega, theta2, T, Rx, Ry, Rz) +
      Rz * compute_pq(std::array<int, 4>
		      ({pqm[0], pqm[1], pqm[2] - 1, pqm[3] + 1}),
		      ints, omega, theta2, T, Rx, Ry, Rz) -
      (pqm[0] - 1) * compute_pq(std::array<int, 4>
				({pqm[0] - 2, pqm[1], pqm[2], pqm[3] + 1}),
				ints, omega, theta2, T, Rx, Ry, Rz) -
      (pqm[1] - 1) * compute_pq(std::array<int, 4>
				({pqm[0], pqm[1] - 2, pqm[2], pqm[3] + 1}),
				ints, omega, theta2, T, Rx, Ry, Rz) -
      (pqm[2] - 1) * compute_pq(std::array<int, 4>
				({pqm[0], pqm[1], pqm[2] - 2, pqm[3] + 1}),
				ints, omega, theta2, T, Rx, Ry, Rz);
    ints[pqm] = res;
    return res;
  }
}

static double compute_apquv(std::array<int, 11> index,
			    std::map<std::array<int, 11>, double> &ints,
			    const std::map<std::array<int, 5>, double> &pquv,
			    std::array<double, 3> c1,
			    std::array<double, 3> c2) {
  if(ints.find(index) != ints.end()) {
    return ints[index];
  } else if(index[0] < 0 || index[1] < 0 || index[2] < 0) {
    return 0;
  } else if(index[0] == 0 && index[1] == 0 && index[2] == 0) {
    std::array<int, 5> ind2 {index[3] + index[6],
      index[4] + index[7],
      index[5] + index[8],
      index[9], index[11]};
    return pquv[ind2] * ((index[7] + index[8] + index[9]) % 2? -1: 1);
  } else {
    std::array<int, 11> ind1 = {index[0] - 1, index[1], index[2],
      index[3] - 1, index[4], index[5],
      index[6], index[7], index[8],
      index[9], index[10]},
      ind2 = {index[0], index[1] - 1, index[2],
      index[3], index[4] - 1, index[5],
      index[6], index[7], index[8],
      index[9], index[10]},
      ind3 = {index[0], index[1], index[2] - 1,
      index[3], index[4], index[5] - 1,
      index[6], index[7], index[8],
      index[9], index[10]},
      ind4 = {index[0] - 1, index[1], index[2],
      index[3], index[4], index[5],
      index[6], index[7], index[8],
      index[9] + 1, index[10] + 1},
      ind5 = {index[0], index[1] - 1, index[2],
      index[3], index[4], index[5],
      index[6], index[7], index[8],
      index[9] + 1, index[10] + 1},
      ind6 = {index[0], index[1], index[2] - 1,
      index[3], index[4], index[5],
      index[6], index[7], index[8],
      index[9] + 1, index[10] + 1},
      ind7 = {index[0] - 1, index[1], index[2],
      index[3] + 1, index[4], index[5],
      index[6], index[7], index[8],
      index[9], index[10] + 1},
      ind8 = {index[0], index[1] - 1, index[2],
      index[3], index[4] + 1, index[5],
      index[6], index[7], index[8],
      index[9], index[10] + 1},
      ind9 = {index[0], index[1], index[2] - 1,
      index[3], index[4], index[5] + 1,
      index[6], index[7], index[8],
      index[9], index[10] + 1};
    double res = index[3] * compute_apquv(ind1, ints, pquv, c1, c2) +
      index[4] * compute_apquv(ind2, ints, pquv, c1, c2) +
      index[5] * compute_apquv(ind3, ints, pquv, c1, c2) -
      (c1[0] - c2[0]) * compute_apquv(ind4, ints, pquv, c1, c2) -
      (c1[1] - c2[1]) * compute_apquv(ind5, ints, pquv, c1, c2) -
      (c1[2] - c2[2]) * compute_apquv(ind6, ints, pquv, c1, c2) +
      compute_apquv(ind7, ints, pquv, c1, c2) +
      compute_apquv(ind8, ints, pquv, c1, c2) +
      compute_apquv(ind9, ints, pquv, c1, c2);
    ints[index] = res;
    return res;
  }
}

static double compute_abcd(std::array<int, 12> index,
			   std::map<std::array<int, 12>, double> &ints,
			   std::array<double, 3> c1,
			   std::array<double, 3> c2,
			   std::array<double, 3> c3,
			   std::array<double, 3> c4) {
  if(ints.find(index) != ints.end()) {
    return ints[index];
  }
  for(int i = 0; i < 12; i++) {
    if(index[i] < 0) {
      return 0;
    }
  }

  if(index[3] == 0 && index[4] == 0 && index[5] == 0 &&
     (index[9] > 0 || index[10] > 0 || index[11] > 0)) {
    std::array<int, 12> ind1 = {
      index[0], index[1], index[2],
      index[3], index[4], index[5],
      index[6] + 1, index[7], index[8],
      index[9] - 1, index[10], index[11]},
      ind2 = {
      index[0], index[1], index[2],
      index[3], index[4], index[5],
      index[6], index[7] + 1, index[8],
      index[9], index[10] - 1, index[11]},
      ind3 = {
      index[0], index[1], index[2],
      index[3], index[4], index[5],
      index[6], index[7], index[8] + 1,
      index[9], index[10], index[11] - 1},
      ind4 = {
      index[0], index[1], index[2],
      index[3], index[4], index[5],
      index[6], index[7], index[8],
      index[9] - 1, index[10], index[11]},
      ind5 = {
      index[0], index[1], index[2],
      index[3], index[4], index[5],
      index[6], index[7], index[8],
      index[9], index[10] - 1, index[11]},
      ind6 = {
      index[0], index[1], index[2],
      index[3], index[4], index[5],
      index[6], index[7], index[8],
      index[9], index[10], index[11] - 1};
    double res = compute_abcd(ind1, ints, c1, c2, c3, c4) +
      compute_abcd(ind2, ints, c1, c2, c3, c4) +
      compute_abcd(ind3, ints, c1, c2, c3, c4) +
      (c3[0] - c4[0]) * compute_abcd(ind4, ints, c1, c2, c3, c4) +
      (c3[1] - c4[1]) * compute_abcd(ind5, ints, c1, c2, c3, c4) +
      (c3[2] - c4[2]) * compute_abcd(ind6, ints, c1, c2, c3, c4);
    ints[index] = res;
    return res;
  } else {
    std::array<int, 12> ind1 = {
      index[0] + 1, index[1], index[2],
      index[3] - 1, index[4], index[5],
      index[6], index[7], index[8],
      index[9], index[10], index[11]},
      ind2 = {
      index[0], index[1] + 1, index[2],
      index[3], index[4] - 1, index[5],
      index[6], index[7], index[8],
      index[9], index[10], index[11]},
      ind3 = {
      index[0], index[1], index[2] + 1,
      index[3], index[4], index[5] - 1,
      index[6], index[7], index[8],
      index[9], index[10], index[11]},
      ind4 = {
      index[0], index[1], index[2],
      index[3] - 1, index[4], index[5],
      index[6], index[7], index[8],
      index[9], index[10], index[11]},
      ind5 = {
      index[0], index[1], index[2],
      index[3], index[4] - 1, index[5],
      index[6], index[7], index[8],
      index[9], index[10], index[11]},
      ind6 = {
      index[0], index[1], index[2],
      index[3], index[4], index[5] - 1,
      index[6], index[7], index[8],
      index[9], index[10], index[11]};
    double res = compute_abcd(ind1, ints, c1, c2, c3, c4) +
      compute_abcd(ind2, ints, c1, c2, c3, c4) +
      compute_abcd(ind3, ints, c1, c2, c3, c4) +
      (c1[0] - c2[0]) * compute_abcd(ind4, ints, c1, c2, c3, c4) +
      (c1[1] - c2[1]) * compute_abcd(ind5, ints, c1, c2, c3, c4) +
      (c1[2] - c2[2]) * compute_abcd(ind6, ints, c1, c2, c3, c4);
    return res;
  }
}

static double integral(const int *pows1, const int *pows2, const int *pows3,
		       const int *pows4, std::array<double, 3> c1,
		       std::array<double, 3> c2,
		       std::array<double, 3> c3,
		       std::array<double, 3> c4,
		       const GaussianOrbital *o1,
		       const GaussianOrbital *o2,
		       const GaussianOrbital *o3,
		       const GaussianOrbital *o4) {
  std::map<std::array<int, 8>, double> e0quv;

  for(int i = 0; i < o3->getnterms(); i++) {
    for(int j = 0; j < o4->getnterms(); j++) {

      std::map<std::array<int, 5>, double> pquv;
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

      for(int k = 0; k < o1->getnterms(); k++) {
	for(int l = 0; l < o2->getnterms(); l++) {
	  double zeta = o1->getalpha(k) + o2->getalpha(l);
	  double Px = (o1->getalpha(k) * c1[0] + o2->getalpha(l) * c2[0]) / zeta,
	    Py = (o1->getalpha(k) * c1[1] + o2->getalpha(l) * c2[1]) / zeta,
	    Pz = (o1->getalpha(k) * c1[2] + o2->getalpha(l) * c2[2]) / zeta;
	  double Kp = M_SQRT2 * std::pow(M_PI, 1.25) / zeta *
	    std::exp(-o1->getalpha(k) * o2->getalpha(l) / eta *
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
	    (std::pow(2 * zeta, o1->getl() + o2->getl()) *
	     std::pow(2 * eta, o3->getl() + o4->getl()) *
	     std::sqrt(zeta + eta));

	  std::map<std::array<int, 4>, double> pqm;

	  // Compute (p|q]_uv
	  for(int px = 0; px <= pows1[0] + pows2[0]; px++) {
	    for(int py = 0; py <= pows1[1] + pows2[1]; py++) {
	      for(int pz = 0; pz <= pows1[2] + pows2[2]; pz++) {
		for(int u = 0; u <= o1->getl() + o2->getl(); u++) {
		  for(int v = 0; v <= o1->getl() + o2->getl(); v++) {
		    std::array<int, 5> index = {px + pows3[0] + pows4[0],
		      py + pows3[1] + pows4[1],
		      pz + pows3[2] + pows4[2],
		      u, v};
		    std::array<int, 4> index2 = {px + pows3[0] + pows4[0],
		      py + pows3[1] + pows4[1],
		      pz + pows3[2] + pows4[2],
		      0};
		    if(pquv.find(index) == pquv.end()) {
		      pquv[index] = 0;
		    }
		    pquv[index] += std::pow(2 * o2->getalpha(l), u) /
		    std::pow(2 * zeta, v) * compute_pq(index2, *&pqm,
						       omega, theta2,
						       T, Rx, Ry, Rz);
		  }
		}
	      }
	    }
	  }
	}
      }
      
      std::map<std::array<int, 11>, double> apquv;

      // Compute (e0|q]
      for(int ex = 0; ex <= pows1[0] + pows2[0]; ex++) {
	for(int ey = 0; ey <= pows1[1] + pows2[1]; ey++) {
	  for(int ez = 0; ez <= pows1[2] + pows2[2]; ez++) {
	    for(int qx = 0; qx <= pows3[0] + pows4[0]; qx++) {
	      for(int qy = 0; qy <= pows3[1] + pows4[1]; qy++) {
		for(int qz = 0; qz <= pows3[2] + pows4[q]; qz++) {
		  for(int u = 0; u <= o3->getl() + o4->getl(); u++) {
		    for(int v = 0; v <= o3->getl() + o4->getl(); v++) {
		      std::array<int, 11> index1 = {ex, ey, ez,
			0, 0, 0,
			qx,
			qy,
			qz,
			0, 0};
		      std::array<int, 8> index2 = {ex, ey, ez,
			qx, qy, qz, u, v};
		      if(e0quv.find(index2) == e0quv.end()) {
			e0quv[index2] = 0;
		      }
		      e0quv[index2] += std::pow(2 * o4->getalpha(j), u) /
			std::pow(2 * eta, v) * compute_apquv(index1, *&apquv,
							     *&pquv, c1, c2);
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

  // Make (e0|f0), then 
      

      



double NumericIntegral::exchange(const GaussianOrbital *o1,
				 const GaussianOrbital *o2,
				 const GaussianOrbital *o3,
				 const GaussianOrbital *o4,
				 std::array<double, 3> c1,
				 std::array<double, 3> c2,
				 std::array<double, 3> c3,
				 std::array<double, 3> c4) const {

}
