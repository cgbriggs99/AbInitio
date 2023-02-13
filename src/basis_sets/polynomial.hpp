
#ifndef __POLYNOMIAL_HPP__
#define __POLYNOMIAL_HPP__

#include <vector>
#include <array>
#include "../util/copyable.hpp"

namespace compchem {
  template<int n>
  class Polynomial : public Copyable<Polynomial<n> > {
  private :
    std::vector<std::array<int, n> > pows;
    std::vector<double> coefs;
  protected :
    void reduce();
  public :
    explicit Polynomial(double scalar);
    Polynomial(const std::vector<std::array<int, n> > &pows, const std::vector<double> &coefs);

    Polynomial(const Polynomial<n> &copy);
    Polynomial();

    virtual ~Polynomial() = default;

    const std::array<int, n> &gettermorder(int index) const;
    double getcoef(int index) const;
    int getsize() const;

    Polynomial<n> &copy() const override;

    Polynomial<n> &operator+=(double rh);
    Polynomial<n> &operator-=(double rh);
    Polynomial<n> &operator*=(double rh);
    Polynomial<n> &operator/=(double rh);
    Polynomial<n> &operator+=(const Polynomial<n> &rh);
    Polynomial<n> &operator-=(const Polynomial<n> &rh);
    Polynomial<n> &operator*=(const Polynomial<n> &rh);
    Polynomial<n> &operator-();
    Polynomial<n> &operator=(double d);
    bool operator==(double d) const;
    bool operator==(const Polynomial<n> &rh) const;
    bool operator!=(double d) const;
    bool operator!=(const Polynomial<n> &rh) const;

    double eval(double x,...) const;

    Polynomial<n> &translate(double x, ...);

    template<int m>
    friend Polynomial<m> &operator+(double lh, const Polynomial<m> &rh);
    template<int m>
    friend Polynomial<m> &operator+(Polynomial<m> lh, double rh);
    template<int m>
    friend Polynomial<m> &operator+(Polynomial<m> lh, const Polynomial<m> &rh);
    template<int m>
    friend Polynomial<m> &operator-(double lh, const Polynomial<m> &rh);
    template<int m>
    friend Polynomial<m> &operator-(Polynomial<m> lh, double rh);
    template<int m>
    friend Polynomial<m> &operator-(Polynomial<m> lh, const Polynomial<m> &rh);
    template<int m>
    friend Polynomial<m> &operator*(double lh, const Polynomial<m> &rh);
    template<int m>
    friend Polynomial<m> &operator*(Polynomial<m> lh, double rh);
    template<int m>
    friend Polynomial<m> &operator*(Polynomial<m> lh, const Polynomial<m> &rh);
    template<int m>
    friend Polynomial<m> &operator/(Polynomial<m> lh, double rh);
    template<int m>
    friend Polynomial<m> &pow(const Polynomial<m> &mant, int expo);
    template<int m>
    friend bool operator==(double lh, const Polynomial<m> &rh);
    template<int m>
    friend bool operator!=(double lh, const Polynomial<m> &rh);
  };
  
  //Polynomial<1> &genlaguerre(int n, double alpha);

  // In cartesian coordinates.
  Polynomial<3> &sphereharm(int l, int ml);
}

#include "polynomial.tcc"


#endif
