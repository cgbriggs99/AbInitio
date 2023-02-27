
#ifndef __POLYNOMIAL_HPP__
#define __POLYNOMIAL_HPP__

#include <vector>
#include <array>
#include <initializer_list>

namespace compchem {
  template<int n>
  class Polynomial {
  private :
    double *coefs;
    int *pows;
    int size;
    
  protected :
    void reduce();
    
  public :
    explicit Polynomial(double scalar);
    Polynomial(const std::vector<std::array<int, n> > &pows, const std::vector<double> &coefs);
    Polynomial(const int *pows, const double *coefs, int size);

    Polynomial(const Polynomial<n> &copy);
    Polynomial();

    virtual ~Polynomial();

    const int *gettermorder(int index) const;
    double getcoef(int index) const;
    int getsize() const;

    void setcoef(int index, double val);

    Polynomial<n> *copy() const;

    Polynomial<n> &operator+=(double rh);
    Polynomial<n> &operator-=(double rh);
    Polynomial<n> &operator*=(double rh);
    Polynomial<n> &operator/=(double rh);
    Polynomial<n> &operator+=(const Polynomial<n> &rh);
    Polynomial<n> &operator-=(const Polynomial<n> &rh);
    Polynomial<n> &operator*=(const Polynomial<n> &rh);
    Polynomial<n> operator-() const;
    bool operator==(double d) const;
    bool operator==(const Polynomial<n> &rh) const;
    bool operator!=(double d) const;
    bool operator!=(const Polynomial<n> &rh) const;

    bool almost(double d, double conv = 1e-6) const;
    bool almost(const Polynomial<n> &rh, double conv = 1e-6) const;

    // Assignment operators.
    Polynomial<n> &operator=(double d);
    Polynomial<n> &operator=(const Polynomial<n> &rh);
    Polynomial<n> &operator=(Polynomial<n> &&rh);

    double eval(double x,...) const;
    double eval(const double *r) const;
    double eval(std::initializer_list<double> r) const;

    Polynomial<n> &translate(double x, ...);

    template<int m>
    friend Polynomial<m> operator+(double lh, const Polynomial<m> &rh);
    template<int m>
    friend Polynomial<m> operator+(Polynomial<m> lh, double rh);
    template<int m>
    friend Polynomial<m> operator+(Polynomial<m> lh, const Polynomial<m> &rh);
    template<int m>
    friend Polynomial<m> operator-(double lh, const Polynomial<m> &rh);
    template<int m>
    friend Polynomial<m> operator-(Polynomial<m> lh, double rh);
    template<int m>
    friend Polynomial<m> operator-(Polynomial<m> lh, const Polynomial<m> &rh);
    template<int m>
    friend Polynomial<m> operator*(double lh, const Polynomial<m> &rh);
    template<int m>
    friend Polynomial<m> operator*(Polynomial<m> lh, double rh);
    template<int m>
    friend Polynomial<m> operator*(Polynomial<m> lh, const Polynomial<m> &rh);
    template<int m>
    friend Polynomial<m> operator/(Polynomial<m> lh, double rh);
    template<int m>
    friend Polynomial<m> pow(const Polynomial<m> &mant, int expo);
    template<int m>
    friend bool operator==(double lh, const Polynomial<m> &rh);
    template<int m>
    friend bool operator!=(double lh, const Polynomial<m> &rh);
  };

  
  
  //Polynomial<1> &genlaguerre(int n, double alpha);

  // In cartesian coordinates.
  Polynomial<3> *sphereharm(int l, int ml);
}

#include "polynomial.template.cpp"


#endif
