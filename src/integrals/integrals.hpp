/*
 * ints.hpp
 *
 *  Created on: Nov 8, 2022
 *      Author: connor
 */

#ifndef GAUSSIAN_HPP_
#define GAUSSIAN_HPP_

#include <array>
#include <vector>
#include <initializer_list>
#include "../util/molecule.hpp"
#include "../util/tei_array.hpp"

namespace compchem {

class IntegralMethod {
public:
  IntegralMethod() = default;
  virtual ~IntegralMethod() = default;

  virtual double overlap(const BasisOrbital *o1, const BasisOrbital *o2,
			 std::array<double, 3> center1,
			 std::array<double, 3> center2) const = 0;
  
  virtual double laplacian(const BasisOrbital *o1, const BasisOrbital *o2,
			   std::array<double, 3> center1,
			   std::array<double, 3> center2) const = 0;
  
  virtual double coulomb(const BasisOrbital *o1, const BasisOrbital *o2,
			 std::array<double, 3> center1,
			 std::array<double, 3> center2,
			 const Atom &atom) const = 0;
  
  virtual double exchange(const BasisOrbital *o1, const BasisOrbital *o2,
			  const BasisOrbital *o3, const BasisOrbital *o4,
			  std::array<double, 3> c1,
			  std::array<double, 3> c2, std::array<double, 3> c3,
			  std::array<double, 3> c4) const = 0;

};

class NumericIntegral : public IntegralMethod {
private :
  int points;
public :
  NumericIntegral(int points) : points(points) {;}
  virtual ~NumericIntegral() = default;

  virtual double overlap(const BasisOrbital *o1, const BasisOrbital *o2,
			 std::array<double, 3> center1,
			 std::array<double, 3> center2) const override;

  virtual double laplacian(const BasisOrbital *o1, const BasisOrbital *o2,
			 std::array<double, 3> center1,
			 std::array<double, 3> center2) const override;

  virtual double coulomb(const BasisOrbital *o1, const BasisOrbital *o2,
			 std::array<double, 3> center1,
			 std::array<double, 3> center2,
			 const Atom &atom) const override;
  
  virtual double exchange(const BasisOrbital *o1, const BasisOrbital *o2,
			  const BasisOrbital *o3, const BasisOrbital *o4,
			  std::array<double, 3> c1,
			  std::array<double, 3> c2, std::array<double, 3> c3,
			  std::array<double, 3> c4) const override;
};
  
  
class AnalyticIntegral {
public :
  AnalyticIntegral() = default;
  virtual ~AnalyticIntegral() = default;

  virtual double overlap(const GaussianOrbital *o1, const GaussianOrbital *o2,
			 std::array<double, 3> center1,
			 std::array<double, 3> center2) const;

  virtual double laplacian(const GaussianOrbital *o1, const GaussianOrbital *o2,
			 std::array<double, 3> center1,
			 std::array<double, 3> center2) const;

  virtual double coulomb(const GaussianOrbital *o1, const GaussianOrbital *o2,
			 std::array<double, 3> center1,
			 std::array<double, 3> center2,
			 const Atom &atom) const;
  
  virtual double exchange(const GaussianOrbital *o1, const GaussianOrbital *o2,
			  const GaussianOrbital *o3, const GaussianOrbital *o4,
			  std::array<double, 3> c1,
			  std::array<double, 3> c2, std::array<double, 3> c3,
			  std::array<double, 3> c4) const;
};  
  
  
template<typename Ints>
class IntegralFactory {
public:
  static void Smatrix(const Molecule *mol, double *out, int *dim);
  static void Tmatrix(const Molecule *mol, double *out, int *dim);
  static void Vmatrix(const Molecule *mol, double *out, int *dim);
  static TEIArray *TEIints(const Molecule *mol);
};

}

#include "integral_factory.template.cpp"

#endif /* GAUSSIAN_HPP_ */
