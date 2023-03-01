 
#include <extramath.h>
#include "../basis_sets/basis_sets.hpp"
#include "../util/atom.hpp"
#include <array>
#include <cmath>
#include "../util/polynomial.hpp"
#include "integrals.hpp"

using namespace compchem;

double AnalyticIntegral::overlap(const GaussianOrbital *o1,
				 const GaussianOrbital *o2,
				 std::array<double, 3> center1,
				 std::array<double, 3> center2) const {

  Polynomial<3> *trans = o2->getharms().copy();
  Polynomial<3> poly = o1->getharms() * trans->translate(center1[0] -center2[0],
							 center1[1] -center2[1],
							 center1[2] -center2[2]);
  delete trans;
  double sum = 0;
  for(int i = 0; i < o1->getnterms(); i++) {
    for(int j = 0; j < o2->getnterms(); j++) {
      // Complete the square.
      double a = o1->getalpha(i) + o2->getalpha(j),
	bx = 2 * o2->getalpha(j) * (center2[0] - center1[0]),
	by = 2 * o2->getalpha(j) * (center2[1] - center1[1]),
	bz = 2 * o2->getalpha(j) * (center2[2] - center1[2]),
	cx = bx * bx / 4 / a,
	cy = by * by / 4 / a,
	cz = bz * bz / 4 / a;
      double coef = o1->getnorm(i) * o2->getnorm(j) *
	std::exp(cx + cy + cz -
		 (center1[0] - center2[0]) *
		 (center1[0] - center2[0]) +
		 (center1[1] - center2[1]) *
		 (center1[1] - center2[1]) +
		 (center1[2] - center2[2]) *
		 (center1[2] - center2[2])) *
	o1->getcoef(i) * o2->getcoef(j);
      // Translate.
      trans = new Polynomial<3>(poly);
      trans->translate(-bx / 2 / a, -by / 2 / a, -bz / 2 / a);
      for(int k = 0; k < trans->getsize(); k++) {
	int pows[3] = trans->gettermorder(k);
	// If a power is odd, then the integral is zero.
	if(pows[0] % 2 || pows[1] % 2 || pows[2] % 2) {
	  continue;
	}
	sum += coef * trans->getcoef(k) * std::tgamma((pows[0] + 1) / 2.0) *
	  std::tgamma((pows[1] + 1) / 2.0) * std::tgamma((pows[2] + 1) / 2.0) /
	  std::pow(a, (pows[0] + pows[1] + pows[2] + 3) / 2.0);
      }
      delete trans;
    }
  }
  return sum;
}
	  
double AnalyticIntegral::laplacian(const GaussianOrbital *o1,
				   const GaussianOrbital *o2,
				   std::array<double, 3> center1,
				   std::array<double, 3> center2) const {
  
  Polynomial<3> *trans = new Polynomial<3>(*o2);
  trans->translate(center1[0] -center2[0],
		   center1[0] -center2[1],
		   center1[0] -center2[2]);
  // Compute the laplacian of the polynomial.
  Polynomial<3> poly = Polynomial<3>(*trans);
  delete trans;
  double one = 1;
  int onepow[5] = {0, 0, 1, 0, 0};
  // Set up coordinate polynomials.
  Polynomial<3> x = Polynomial<3>(onepow + 2, &one, 1),
    y = Polynomial<3>(onepow + 1, &one, 1),
    z = Polynomial<3>(onepow, &one, 1);
  
  double sum = 0;
  for(int i = 0; i < o1->getnterms(); i++) {
    Polynomial<3> lap = Polynomial<3>(0);    // This function is symmetric.

    // Compute the laplacian for the term.
    for(int j = 0; j < o1->getharms().getsize(); j++) {
      int pows1[3] = o1->getharms().gettermorder(j);
      lap += o1->getharms().getcoef(j) *
	((pows1[0] >= 2? (pows1[0] * (pows1[0] - 1) *
			  compchem::pow(x, pows1[0] - 2) *
			  compchem::pow(y, pows1[1]) *
			  compchem::pow(z, pows1[2])) : 0) +
	 (pows1[1] >= 2? (pows1[1] * (pows1[1] - 1) *
			  compchem::pow(x, pows1[0]) *
			  compchem::pow(y, pows[1] - 2) *
			  compchem::pow(z, pows[2])) : 0) +
	 (pows1[2] >= 2? (pows1[2] * (pows1[2] - 1) *
			  compchem::pow(x, pows1[0]) *
			  compchem::pow(y, pows1[1]) *
			  compchem::pow(z, pows1[2] - 2)) : 0) +
	 4 * o1->getalpha(i) * ((pows1[0] >= 1?
				 (pows1[0] * compchem::pow(x, pows1[0]) *
				  compchem::pow(y, pows1[1]) *
				  compchem::pow(z, pows1[2])): 0) +
				(pows1[1] >= 1?
				 (pows1[1] * compchem::pow(x, pows1[0]) *
				  compchem::pow(y, pows1[1]) *
				  compchem::pow(z, pows1[2])): 0) +
				(pows1[2] >= 1?
				 (pows1[2] * compchem::pow(x, pows1[0]) *
				  compchem::pow(y, pows1[1]) *
				  compchem::pow(z, pows1[2])): 0)) +
	 compchem::pow(x, pows1[0]) *
	 compchem::pow(y, pows1[1]) *
	 compchem::pow(z, pows1[2]) *
	 (4 * (x * x + y * y + z * z) o1->getalpha(i) * o1->getalpha(i) +
	  2 * o1->getalpha(i)));
    }
							  
    for(int j = 0; j < o2->getnterms(); j++) {
      // Complete the square.
      double a = o1->getalpha(i) + o2->getalpha(j),
	bx = 2 * o2->getalpha(j) * (center2[0] - center1[0]),
	by = 2 * o2->getalpha(j) * (center2[1] - center1[1]),
	bz = 2 * o2->getalpha(j) * (center2[2] - center1[2]),
	cx = bx * bx / 4 / a,
	cy = by * by / 4 / a,
	cz = bz * bz / 4 / a;
      double coef = o1->getnorm(i) * o2->getnorm(j) *
	std::exp(cx + cy + cz -
		 (center1[0] - center2[0]) *
		 (center1[0] - center2[0]) +
		 (center1[1] - center2[1]) *
		 (center1[1] - center2[1]) +
		 (center1[2] - center2[2]) *
		 (center1[2] - center2[2])) *
	o1->getcoef(i) * o2->getcoef(j);
      // Translate.
      trans = new Polynomial<3>(lap * poly);
      trans->translate(-bx / 2 / a, -by / 2 / a, -bz / 2 / a);
      for(int k = 0; k < trans->getsize(); k++) {
	int pows2[3] = trans->gettermorder(k);
	// If a power is odd, then the integral is zero.
	if(pows2[0] % 2 || pows2[1] % 2 || pows2[2] % 2) {
	  continue;
	}
	sum += coef * trans->getcoef(k) * std::tgamma((pows2[0] + 1) / 2.0) *
	  std::tgamma((pows2[1] + 1) / 2.0) * std::tgamma((pows2[2] + 1) / 2.0) /
	  std::pow(a, (pows2[0] + pows2[1] + pows2[2] + 3) / 2.0);
      }
      delete trans;
    }
  }
  return sum;
}
