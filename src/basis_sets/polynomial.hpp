
#ifndef __POLYNOMIAL_HPP__
#define __POLYNOMIAL_HPP__

namespace compchem {
  template<int n>
  class Polynomial : Copyable {
  private :
    std::vector<std::array<int, n> > pows;
    std::vector<double> coefs;
  public:
    Polynomial(double scalar);
    Polynomial(std::vector<std::array<int, n> > pows, std::vector<double> coefs);
    Polynomial(std::initializer_list<int> pows, double coef);

    Polynomial(const Polynomial &copy);
    Polynomial();

    virtual ~Polynomial() = default;

    const std::array<int, n> &gettermorder(int index) const;
    double getcoef(int index) const;
    int getsize() const;

    Polynomial<n> &copy() override;

    Polynomial<n> &operator+=(const Polynomial<n> &rh);
    Polynomial<n> &operator-=(const Polynomial<n> &rh);
    Polynomial<n> &operator*=(const Polynomial<n> &rh);
    Polynomial<n> &operator+=(double rh);
    Polynomial<n> &operator-=(double rh);
    Polynomial<n> &operator*=(double rh);
    Polynomial<n> &operator/=(double rh);
    Polynomial<n> &operator-();

    double operator()(double x, ...);

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
  };

  template<int m>
  Polynomial<m> &operator+(double lh, const Polynomial<m> &rh);
  template<int m>
  Polynomial<m> &operator+(Polynomial<m> lh, double rh);
  template<int m>
  Polynomial<m> &operator+(Polynomial<m> lh, const Polynomial<m> &rh);
  template<int m>
  Polynomial<m> &operator-(double lh, const Polynomial<m> &rh);
  template<int m>
  Polynomial<m> &operator-(Polynomial<m> lh, double rh);
  template<int m>
  Polynomial<m> &operator-(Polynomial<m> lh, const Polynomial<m> &rh);
  template<int m>
  Polynomial<m> &operator*(double lh, const Polynomial<m> &rh);
  template<int m>
  Polynomial<m> &operator*(Polynomial<m> lh, double rh);
  template<int m>
  Polynomial<m> &operator*(Polynomial<m> lh, const Polynomial<m> &rh);
  template<int m>
  Polynomial<m> &operator/(Polynomial<m> lh, double rh);

  Polynomial<1> &genlaguerre(int n, double alpha);

  // In cartesian coordinates.
  Polynomial<3> &sphereharm(int l, int ml);
}


#endif
