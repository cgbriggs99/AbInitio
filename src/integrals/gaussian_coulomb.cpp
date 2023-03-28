 
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

double AnalyticIntegral::os_coul(const std::array<int, 7> &index,
				 std::map<std::array<int, 7>, double> &ints,
				 const std::array<double, 3> &c1,
				 const std::array<double, 3> &c2,
				 const std::array<double, 3> &ca,
				 double Rx, double Ry, double Rz,
				 double Px, double Py, double Pz,
				 double zeta) const {
  if(ints.count(index) != 0) {
    return ints[index];
  }

  for(int i = 0; i < 7; i++) {
    if(index[i] < 0) {
      return 0;
    }
  }

  if(index[0] == 0) {
    if(index[1] == 0) {
      if(index[2] == 0) {
	if(index[3] == 0) {
	  if(index[4] == 0) {
	    if(index[5] == 0) {
	      double res = boys_square(index[6], zeta *
					(Rx * Rx + Ry * Ry + Rz * Rz));
	      ints[index] = res;
	      return res;
	    } else {
	      std::array<int, 7> ind1 = {
		index[0], index[1], index[2],
		index[3], index[4], index[5] - 1,
		index[6]},
		ind2 = {
		index[0], index[1], index[2],
		index[3], index[4], index[5] - 2,
		index[6]},
		ind3 = {
		index[0], index[1], index[2],
		index[3], index[4], index[5] - 1,
		index[6] + 1},
		ind4 = {
		index[0], index[1], index[2],
		index[3], index[4], index[5] - 2,
		index[6] + 1};
	      double res = (Pz - c2[2]) * os_coul(ind1, ints, c1, c2, ca,
						  Px, Py, Pz,
						  Rx, Ry, Rz,
						  zeta) +
		(index[5] - 1) * os_coul(ind2, ints, c1, c2, ca,
					 Px, Py, Pz,
					 Rx, Ry, Rz,
					 zeta) / (2 * zeta) -
		(Pz - ca[2]) * os_coul(ind3, ints, c1, c2, ca,
				       Px, Py, Pz,
				       Rx, Ry, Rz,
				       zeta) -
		(index[5] - 1) * os_coul(ind4, ints, c1, c2, ca,
					 Px, Py, Pz,
					 Rx, Ry, Rz,
					 zeta) / (2 * zeta);
	      ints[index] = res;
	      return res;
	    }
	  } else {
	    std::array<int, 7> ind1 = {
	      index[0], index[1], index[2],
	      index[3], index[4] - 1, index[5],
	      index[6]},
	      ind2 = {
		index[0], index[1], index[2],
		index[3], index[4] - 2, index[5],
		index[6]},
	      ind3 = {
		index[0], index[1], index[2],
		index[3], index[4] - 1, index[5],
		index[6] + 1},
	      ind4 = {
		index[0], index[1], index[2],
		index[3], index[4] - 2, index[5],
		index[6] + 1};
	    double res = (Py - c2[1]) * os_coul(ind1, ints, c1, c2, ca,
						Px, Py, Pz,
						Rx, Ry, Rz,
						zeta) +
	      (index[4] - 1) * os_coul(ind2, ints, c1, c2, ca,
				       Px, Py, Pz,
				       Rx, Ry, Rz,
				       zeta) / (2 * zeta) -
	      (Py - ca[1]) * os_coul(ind3, ints, c1, c2, ca,
				     Px, Py, Pz,
				     Rx, Ry, Rz,
				     zeta) -
	      (index[4] - 1) * os_coul(ind4, ints, c1, c2, ca,
				       Px, Py, Pz,
				       Rx, Ry, Rz,
				       zeta) / (2 * zeta);
	    ints[index] = res;
	    return res;
	  }
	} else {
	  std::array<int, 7> ind1 = {
	    index[0], index[1], index[2],
	    index[3] - 1, index[4], index[5],
	    index[6]},
	    ind2 = {
	      index[0], index[1], index[2],
	      index[3] - 2, index[4], index[5],
	      index[6]},
	    ind3 = {
	      index[0], index[1], index[2],
	      index[3] - 1, index[4], index[5],
	      index[6] + 1},
	    ind4 = {
	      index[0], index[1], index[2],
	      index[3] - 2, index[4], index[5],
	      index[6] + 1};
	  double res = (Px - c2[0]) * os_coul(ind1, ints, c1, c2, ca,
					      Px, Py, Pz,
					      Rx, Ry, Rz,
					      zeta) +
	    (index[3] - 1) * os_coul(ind2, ints, c1, c2, ca,
				     Px, Py, Pz,
				     Rx, Ry, Rz,
				     zeta) / (2 * zeta) -
	    (Px - ca[0]) * os_coul(ind3, ints, c1, c2, ca,
				   Px, Py, Pz,
				   Rx, Ry, Rz,
				   zeta) -
	    (index[3] - 1) * os_coul(ind4, ints, c1, c2, ca,
				     Px, Py, Pz,
				     Rx, Ry, Rz,
				     zeta) / (2 * zeta);
	  ints[index] = res;
	  return res;
	}
      } else {
	std::array<int, 7> ind1 = {
	  index[0], index[1], index[2] - 1,
	  index[3], index[4], index[5],
	  index[6]},
	  ind2 = {
	    index[0], index[1], index[2] - 2,
	    index[3], index[4], index[5],
	    index[6]},
	  ind3 = {
	    index[0], index[1], index[2] - 1,
	    index[3], index[4], index[5] - 1,
	    index[6]},
	  ind4 = {
	    index[0], index[1], index[2] - 1,
	    index[3], index[4], index[5],
	    index[6] + 1},
	  ind5 = {
	    index[0], index[1], index[2] - 2,
	    index[3], index[4], index[5],
	    index[6] + 1},
	  ind6 = {
	    index[0], index[1], index[2] - 1,
	    index[3], index[4], index[5] - 1,
	    index[6] + 1};
	double res = (Pz - c1[2]) * os_coul(ind1, ints, c1, c2, ca,
					    Px, Py, Pz,
					    Rx, Ry, Rz,
					    zeta) +
	  (index[2] - 1) * os_coul(ind2, ints, c1, c2, ca,
				   Px, Py, Pz,
				   Rx, Ry, Rz,
				   zeta) / (2 * zeta) + 
	  index[5] * os_coul(ind3, ints, c1, c2, ca,
			     Px, Py, Pz,
			     Rx, Ry, Rz,
			     zeta) / (2 * zeta) -
	  (Pz - ca[2]) * os_coul(ind4, ints, c1, c2, ca,
				 Px, Py, Pz,
				 Rx, Ry, Rz,
				 zeta) -
	  (index[2] - 1) * os_coul(ind5, ints, c1, c2, ca,
				   Px, Py, Pz,
				   Rx, Ry, Rz,
				   zeta) / (2 * zeta) -
	  index[5] * os_coul(ind6, ints, c1, c2, ca,
			     Px, Py, Pz,
			     Rx, Ry, Rz,
			     zeta) / (2 * zeta);
	ints[index] = res;
	return res;
      }
    } else {
      std::array<int, 7> ind1 = {
	index[0], index[1] - 1, index[2],
	index[3], index[4], index[5],
	index[6]},
	ind2 = {
	  index[0], index[1] - 2, index[2],
	  index[3], index[4], index[5],
	  index[6]},
	ind3 = {
	  index[0], index[1] - 1, index[2],
	  index[3], index[4] - 1, index[5],
	  index[6]},
	ind4 = {
	  index[0], index[1] - 1, index[2],
	  index[3], index[4], index[5],
	  index[6] + 1},
	ind5 = {
	  index[0], index[1] - 2, index[2],
	  index[3], index[4], index[5],
	  index[6] + 1},
	ind6 = {
	  index[0], index[1] - 1, index[2],
	  index[3], index[4] - 1, index[5],
	  index[6] + 1};
      double res = (Py - c1[1]) * os_coul(ind1, ints, c1, c2, ca,
					  Px, Py, Pz,
					  Rx, Ry, Rz,
					  zeta) +
	(index[1] - 1) * os_coul(ind2, ints, c1, c2, ca,
				 Px, Py, Pz,
				 Rx, Ry, Rz,
				 zeta) / (2 * zeta) + 
	index[4] * os_coul(ind3, ints, c1, c2, ca,
			   Px, Py, Pz,
			   Rx, Ry, Rz,
			   zeta) / (2 * zeta) -
	(Py - ca[1]) * os_coul(ind4, ints, c1, c2, ca,
			       Px, Py, Pz,
			       Rx, Ry, Rz,
			       zeta) -
	(index[1] - 1) * os_coul(ind5, ints, c1, c2, ca,
				 Px, Py, Pz,
				 Rx, Ry, Rz,
				 zeta) / (2 * zeta) -
	index[4] * os_coul(ind6, ints, c1, c2, ca,
			   Px, Py, Pz,
			   Rx, Ry, Rz,
			   zeta) / (2 * zeta);
      ints[index] = res;
      return res;
    }
  } else {
    std::array<int, 7> ind1 = {
      index[0] - 1, index[1], index[2],
      index[3], index[4], index[5],
      index[6]},
      ind2 = {
	index[0] - 2, index[1], index[2],
	index[3], index[4], index[5],
	index[6]},
      ind3 = {
	index[0] - 1, index[1], index[2],
	index[3] - 1, index[4], index[5],
	index[6]},
      ind4 = {
	index[0] - 1, index[1], index[2],
	index[3], index[4], index[5],
	index[6] + 1},
      ind5 = {
	index[0] - 2, index[1], index[2],
	index[3], index[4], index[5],
	index[6] + 1},
      ind6 = {
	index[0] - 1, index[1], index[2],
	index[3] - 1, index[4], index[5],
	index[6] + 1};
    double res = (Px - c1[0]) * os_coul(ind1, ints, c1, c2, ca,
					Px, Py, Pz,
					Rx, Ry, Rz,
					zeta) +
      (index[0] - 1) * os_coul(ind2, ints, c1, c2, ca,
			       Px, Py, Pz,
			       Rx, Ry, Rz,
			       zeta) / (2 * zeta) + 
      index[3] * os_coul(ind3, ints, c1, c2, ca,
			 Px, Py, Pz,
			 Rx, Ry, Rz,
			 zeta) / (2 * zeta) -
      (Px - ca[0]) * os_coul(ind4, ints, c1, c2, ca,
			     Px, Py, Pz,
			     Rx, Ry, Rz,
			     zeta) -
      (index[0] - 1) * os_coul(ind5, ints, c1, c2, ca,
			       Px, Py, Pz,
			       Rx, Ry, Rz,
			       zeta) / (2 * zeta) -
      index[3] * os_coul(ind6, ints, c1, c2, ca,
			 Px, Py, Pz,
			 Rx, Ry, Rz,
			 zeta) / (2 * zeta);
    ints[index] = res;
    return res;
  }
}
    

double AnalyticIntegral::coul_integral(const int *pows1,
				       const int *pows2,
				       const std::array<double, 3> &c1,
				       const std::array<double, 3> &c2,
				       const GaussianOrbital *o1,
				       const GaussianOrbital *o2,
				       const Atom &atom) const {

  double sum = 0;
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
	//-
	omega = o1->getcoef(i) * o2->getcoef(j) *
	o1->getnorm(i) * o2->getnorm(j) * K;
      std::map<std::array<int, 7>, double> os_ints;
      std::array<int, 7> index = {
	pows1[0], pows1[1], pows1[2],
	pows2[0], pows2[1], pows2[2],
	0};
      sum += omega * os_coul(index, os_ints, c1, c2, atom.getpos(),
			     Rx, Ry, Rz, Px, Py, Pz, zeta);
    }
  }
  return sum;
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
