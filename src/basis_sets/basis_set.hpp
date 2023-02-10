 
#ifndef __BASIS_SET_HPP__
#define __BASIS_SET_HPP__

#include "../util/copyable.hpp"
#include "electron_configs.hpp"
#include <vector>
#include <initializer_list>

namespace compchem {

  class BasisOrbital : Copyable {
  public :
    virtual ~BasisOrbital() = 0;

    virtual double eval(double x, double y, double z) const = 0;
  };

  class SlaterOrbital : BasisOrbital {
  private :
    double Zeff;
    double neff;
    int l, ml;
  public :
    SlaterOrbital(double Zeff, double neff, int l, int ml);
    virtual ~SlaterOrbital() = default;
    SlaterOrbital(const SlaterOrbital &copy);

    double eval(double x, double y, double z) const override;

    SlaterOrbital &copy() override;
  };

  double slater_rule(int n, int l, const compchem::GSConfig &conf);

  class GaussianOrbital : BasisOrbital {
  private :
    std::vector<double> coefs, alphas;
    int l, ml;
  public :
    GaussianOrbital(int l, int ml, const std::vector<double> &coefs, const std::vector<double> &alphas);
    GaussianOrbital(const GaussianOrbital &copy);
    virtual ~GaussianOrbital() = default;

    const std::vector<double> getcoefs() const;
    const std::vector<double> getalphas() const;
    double getcoef(int index) const;
    double getalpha(int index) const;
    int getl() const;
    int getml() const;

    double eval(double x, double y, double z) const override;

    GaussianOrbital &copy() override;
  };

  
  
}


#endif
