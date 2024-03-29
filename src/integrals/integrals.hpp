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
#include "../opts/options.hpp"

namespace compchem {

class IntegralMethod {
protected:
  OptionList &opts;
public:
  IntegralMethod(const OptionList &opts) : opts(*new OptionList(opts)) {;}
  IntegralMethod() : opts(*new OptionList(GlobalOptions::getsingleton())) {;}
  virtual ~IntegralMethod();

  virtual double overlap(const BasisOrbital *o1, const BasisOrbital *o2,
			 std::array<double, 3> center1,
			 std::array<double, 3> center2) const = 0;
  
  virtual double kinetic(const BasisOrbital *o1, const BasisOrbital *o2,
			   std::array<double, 3> center1,
			   std::array<double, 3> center2) const = 0;
  
  virtual double attraction(const BasisOrbital *o1, const BasisOrbital *o2,
			 std::array<double, 3> center1,
			 std::array<double, 3> center2,
			 const Atom &atom) const = 0;
  
  virtual double repulsion(const BasisOrbital *o1, const BasisOrbital *o2,
			  const BasisOrbital *o3, const BasisOrbital *o4,
			  std::array<double, 3> c1,
			  std::array<double, 3> c2, std::array<double, 3> c3,
			  std::array<double, 3> c4) const = 0;

};

class NumericIntegral : public IntegralMethod {
private :
  int points;
public :
  NumericIntegral(int points, const OptionList &opts) : points(points),
						 IntegralMethod(opts) {;}
  NumericIntegral(int points) : points(points),
				IntegralMethod()
  {;}
  virtual ~NumericIntegral() = default;

  virtual double overlap(const BasisOrbital *o1, const BasisOrbital *o2,
			 std::array<double, 3> center1,
			 std::array<double, 3> center2) const override;

  virtual double kinetic(const BasisOrbital *o1, const BasisOrbital *o2,
			 std::array<double, 3> center1,
			 std::array<double, 3> center2) const override;

  virtual double attraction(const BasisOrbital *o1, const BasisOrbital *o2,
			 std::array<double, 3> center1,
			 std::array<double, 3> center2,
			 const Atom &atom) const override;
  
  virtual double repulsion(const BasisOrbital *o1, const BasisOrbital *o2,
			  const BasisOrbital *o3, const BasisOrbital *o4,
			  std::array<double, 3> c1,
			  std::array<double, 3> c2, std::array<double, 3> c3,
			  std::array<double, 3> c4) const override;
};
  
  
class AnalyticIntegral {
protected:
  OptionList &opts;

  double os_attr(const std::array<int, 7> &index,
		 std::map<std::array<int, 7>, double> &ints,
		 const std::array<double, 3> &c1,
		 const std::array<double, 3> &c2,
		 const std::array<double, 3> &ca,
		 double Rx, double Ry, double Rz,
		 double Px, double Py, double Pz,
		 double zeta) const;

  double attr_integral(const int *pows1,
		       const int *pows2,
		       const std::array<double, 3> &c1,
		       const std::array<double, 3> &c2,
		       const GaussianOrbital *o1,
		       const GaussianOrbital *o2,
		       const Atom &atom) const;

  double rep_integral(const int *pows1, const int *pows2, const int *pows3,
		      const int *pows4,
		      const std::array<double, 3> &c1,
		      const std::array<double, 3> &c2,
		      const std::array<double, 3> &c3,
		      const std::array<double, 3> &c4,
		      const GaussianOrbital *o1,
		      const GaussianOrbital *o2,
		      const GaussianOrbital *o3,
		      const GaussianOrbital *o4) const;
  
public :
  AnalyticIntegral(const OptionList &opts) : opts(*new OptionList(opts)) {;}
  AnalyticIntegral() : opts(*new OptionList(GlobalOptions::getsingleton())) {;}
  virtual ~AnalyticIntegral();

  virtual double overlap(const GaussianOrbital *o1, const GaussianOrbital *o2,
			 std::array<double, 3> center1,
			 std::array<double, 3> center2) const;

  virtual double kinetic(const GaussianOrbital *o1, const GaussianOrbital *o2,
			 std::array<double, 3> center1,
			 std::array<double, 3> center2) const;

  virtual double attraction(const GaussianOrbital *o1, const GaussianOrbital *o2,
			    std::array<double, 3> center1,
			    std::array<double, 3> center2,
			    const Atom &atom) const;
  
  virtual double repulsion(const GaussianOrbital *o1, const GaussianOrbital *o2,
			   const GaussianOrbital *o3, const GaussianOrbital *o4,
			   std::array<double, 3> c1,
			   std::array<double, 3> c2, std::array<double, 3> c3,
			   std::array<double, 3> c4) const;

  // Gaussian kernel for the Boys integral. For coulomb.
  double boys_square(int j, double T) const;
  // Exponential kernel for the Boys interal. For exchange.
  double boys_linear(int j, double T) const;
};
  
  
template<typename Ints>
class IntegralFactory {
protected:
  OptionList &opts;
  static void s_routine(const std::vector<const compchem::GaussianOrbital *> *orbs,
			const std::vector<std::array<double, 3> *> *centers,
			double *out,
			int dim, int thread_num, int threads,
			compchem::OptionList &opts);
  static void t_routine(const std::vector<const compchem::GaussianOrbital *> *orbs,
			const std::vector<std::array<double, 3> *> *centers,
			double *out,
			int dim, int thread_num, int threads,
			compchem::OptionList &opts);
  static void v_routine(const std::vector<const compchem::GaussianOrbital *> *orbs,
			const std::vector<std::array<double, 3> *> *centers,
			double *out,
			int dim, int thread_num, int threads,
			compchem::OptionList &opts,
			const compchem::Molecule *mol);

  static void tei_routine(const std::vector<const compchem::GaussianOrbital *> *orbs,
			const std::vector<std::array<double, 3> *> *centers,
			TEIArray *out,
			int dim, int thread_num, int threads,
			compchem::OptionList &opts);
public:
  IntegralFactory() : opts(*new OptionList(GlobalOptions::getsingleton())) {;}
  IntegralFactory(const OptionList &opts) : opts(*new OptionList(opts)) {;}
  virtual ~IntegralFactory();
  
  void Smatrix(const Molecule *mol, double *out, int *dim);
  void Tmatrix(const Molecule *mol, double *out, int *dim);
  void Vmatrix(const Molecule *mol, double *out, int *dim);
  TEIArray *TEIints(const Molecule *mol);
};

}

#include "integral_factory.template.cpp"

#endif /* GAUSSIAN_HPP_ */
