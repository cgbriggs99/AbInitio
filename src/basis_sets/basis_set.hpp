 
#ifndef __BASIS_SET_HPP__
#define __BASIS_SET_HPP__

#include "../util/electron_configs.hpp"
#include "../util/polynomial.hpp"
#include <vector>
#include <initializer_list>
#include <exception>

namespace compchem {

  class BasisOrbital {
  public :
    virtual ~BasisOrbital() = default;

    virtual double eval(double x, double y, double z) const = 0;
    virtual BasisOrbital *copy() const = 0;
  };

  class SlaterOrbital : public BasisOrbital {
  private :
    double Zeff;
    int n, l, ml;
    Polynomial<3> *harms;
  public :
    SlaterOrbital(double Zeff, int n, int l, int ml);
    virtual ~SlaterOrbital();
    SlaterOrbital(const SlaterOrbital &copy);

    double eval(double x, double y, double z) const override;

    double getZeff() const;
    int getn() const;
    int getl() const;
    int getml() const;
    const Polynomial<3> &getharms() const;

    BasisOrbital *copy() const override;
  };

  double slater_rule(int n, int l, const compchem::GSConfig &conf);

  class GaussianOrbital : public BasisOrbital {
  private :
    std::vector<double> coefs, alphas;
    int l, ml;
    Polynomial<3> *harms;
  public :
    GaussianOrbital(int l, int ml, const std::vector<double> &coefs, const std::vector<double> &alphas);
    GaussianOrbital(const GaussianOrbital &copy);
    virtual ~GaussianOrbital();

    const std::vector<double> &getcoefs() const;
    const std::vector<double> &getalphas() const;
    double getcoef(int index) const;
    double getalpha(int index) const;
    int getnterms() const;
    int getl() const;
    int getml() const;
    const Polynomial<3> &getharms() const;

    double eval(double x, double y, double z) const override;

    BasisOrbital *copy() const override;
  };

  
  
}


#endif
