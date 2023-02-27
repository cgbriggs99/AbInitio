
#include <cmath>
#include "polynomial.hpp"

using namespace compchem;

Polynomial<3> *compchem::sphereharm(int l, int ml) {
  std::vector<std::array<int, 3> > r2pow = {{2, 0, 0}, {0, 2, 0}, {0, 0, 2}},
    xpow = {{1, 0, 0}}, ypow = {{0, 1, 0}}, zpow {{0, 0, 1}};
  std::vector<double> r2cos = {1, 1, 1}, xcos = {1};
  Polynomial<3> *pi = new Polynomial<3>(), *r2 = new Polynomial<3>(r2pow, r2cos);
  Polynomial<3> *x = new Polynomial<3>(xpow, xcos),
    *y = new Polynomial<3>(ypow, xcos),
    *z = new Polynomial<3>(zpow, xcos);

  if(ml > 0) {
    Polynomial<3> *a = new Polynomial<3>();
    double norm = sqrt((2 * l + 1) / (2 * M_PI));
    double picoef = exp(0.5 * (-lgamma(l - ml + 1) - lgamma(l + ml + 1)) -
		       lgamma(l + 1) + lgamma(2 * l + 1)) * std::pow(2, -l);
    double acoef = 1;

    for(int p = 0; p <= ml; p++) {
      if((ml - p) % 4 == 0) {
	*a += acoef * compchem::pow(*x, p) * compchem::pow(*y, ml - p);
      } else if((ml - p) % 4 == 2) {
	*a -= acoef * compchem::pow(*x, p) * compchem::pow(*y, ml - p);
      }
      acoef *= (ml - p);
      acoef /= p + 1;
    }

    for(int k = 0; k <= (l - ml) / 2; k++) {
      *pi += picoef * compchem::pow(*r2, 2 * k) * compchem::pow(*z, l - 2 * k - ml);
      picoef *= -(l - k) * (l - 2 * k - ml) * (l - 2 * k - 1 - ml);
      picoef /= (k + 1) * (2 * l - 2 * k) * (2 * l - 2 * k - 1);
    }

    *pi *= norm * *a;
    delete a;
    delete r2;
    delete x;
    delete y;
    delete z;
    return pi;
  } else if(ml < 0) {
    Polynomial<3> *b = new Polynomial<3>();
    double norm = sqrt((2 * l + 1) / (2 * M_PI));
    double picoef = exp(0.5 * (-lgamma(l + ml + 1) - lgamma(l - ml + 1)) -
		       lgamma(l + 1) + lgamma(2 * l + 1)) * std::pow(2, -l);
    double bcoef = 1;

    for(int p = 0; p <= -ml; p++) {
      if((-ml - p) % 4 == 1) {
	*b += bcoef * compchem::pow(*x, p) * compchem::pow(*y, -ml - p);
      } else if((-ml - p) % 4 == 3) {
	*b -= bcoef * compchem::pow(*x, p) * compchem::pow(*y, -ml - p);
      }
      bcoef *= (ml - p);
      bcoef /= p + 1;
    }

    for(int k = 0; k <= (l + ml) / 2; k++) {
      *pi += picoef * compchem::pow(*r2, 2 * k) *
	compchem::pow(*z, l - 2 * k + ml);
      picoef *= -(l - k) * (l - 2 * k + ml) * (l - 2 * k - 1 + ml);
      picoef /= (k + 1) * (2 * l - 2 * k) * (2 * l - 2 * k - 1);
    }

    *pi *= norm * *b;
    delete b;
    delete r2;
    delete x;
    delete y;
    delete z;
    return pi;
  } else {
    double norm = sqrt((2 * l + 1) / (4 * M_PI));
    double picoef = exp(lgamma(2 * l + 1) - 2 * lgamma(l + 1)) *
      std::pow(2, -l);

    for(int k = 0; k <= l / 2; k++) {
      *pi += picoef * compchem::pow(*r2, 2 * k) * compchem::pow(*z, l - 2 * k);
      picoef *= -(l - k) * (l - 2 * k) * (l - 2 * k - 1);
      picoef /= (k + 1) * (2 * l - 2 * k) * (2 * l - 2 * k - 1);
    }

    *pi *= norm;
    delete r2;
    delete x;
    delete y;
    delete z;
    return pi;
  }
}
