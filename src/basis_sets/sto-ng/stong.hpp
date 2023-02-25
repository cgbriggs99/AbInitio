/*
 * sto3g.hpp
 *
 *  Created on: Feb 7, 2023
 *      Author: connor
 */

#ifndef STONG_HPP_
#define STONG_HPP_

#include "../basis_set.hpp"
#include <vector>
#include <cmath>

namespace compchem {
  
class STO_nGComputer {
private :
  int n;
  
  void gradientdescent(std::vector<double> &coefs,
		       std::vector<double> &alphas,
		       double zeff, int N, int l) const;
  
  void gradientdescentptype(std::vector<double> &coefs,
			    const std::vector<double> &alphas,
			    double zeff, int N, int l) const;
  
public:
  ~STO_nGComputer() = default;
  STO_nGComputer(int n) :
    n(n) {;}
  
  int getnGauss() const;
  
  std::vector<GaussianOrbital> *compute(int Z) const;

  double gfunc(int l, double alpha1, double alpha2) const;
  double gderiva1(int l, double alpha1, double alpha2) const;
  double sfunc(int N, int l, double alpha, double zeta) const;
  double sderiv(int N, int l, double alpha, double zeta) const;
  
  double loss(const std::vector<double> &coefs,
	      const std::vector<double> &alphas,
	      double zeff, int N, int l, double lambda) const;
  
  std::vector<double> gradient(const std::vector<double> &coefs,
			       const std::vector<double> &alphas,
			       double zeff, int N, int l, double lambda) const;
  
  std::vector<double> gradientptype(const std::vector<double> &coefs,
				    const std::vector<double> &alphas,
				    double zeff, int N, int l, double lambda) const;
  
};


}

#endif /* STO3G_HPP_ */
