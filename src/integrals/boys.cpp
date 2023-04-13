#include "integrals.hpp"
#include <cmath>
#include "extramath.h"
#include <map>

#define TOL 1e-6

static double integrate_square(int j, double T, int points) {
  double step = 1.0 / (points - 1);
  double sum = 0;
  for(int i = 1; i < points - 1; i++) {
    sum += step * std::pow(i * step, 2 * j) *
      std::exp(-T * (i * i * step * step));
  }
  sum += step * std::exp(-T) / 2;
  return sum;
}

static double extrapolant_square(int i, int k, int j, double T,
			  std::map<std::pair<int, int>, double> &ints) {
  if(ints.count(std::pair<int, int>(i, k)) != 0) {
    return ints.at(std::pair<int, int>(i, k));
  }
  if(i == 0) {
    double res = integrate_square(j, T, 2 << k);
    ints[std::pair<int, int>(i, k)] = res;
    return res;
  } else {
    double res = extrapolant_square(i - 1, k + 1, j, T, ints) -
      extrapolant_square(i - 1, k, j, T, ints) / (std::pow(2, k) - 1);
    ints[std::pair<int, int>(i, k)] = res;
    return res;
  }
}
  

double compchem::AnalyticIntegral::boys_square(int j, double T) const {
  if(j < 0) {
    return 0;
  }
  // If T = 0, this is just a power function.
  if(T == 0) {
    return 1.0 / (2.0 * j + 1.0);
  }
  // If j = 0, this integral is just an error function.
  if(j == 0) {
    return std::erf(std::sqrt(T)) * std::sqrt(M_PI / T) / 2;
  }

  if(this->opts.getbooloption("analytic boys")) {
    if(T < 1) {
      double sum1 = 0, sum2 = 1, term = 1;
      for(int i = 0; i < 30 && 
	    sum2 != sum1; i++) {
	sum2 = sum1;
	sum1 += term / (2 * j + 2 * i + 1);
        term *= -T / (i + 1);
      }
      return sum1;
    } else {
      double sum = 0, term = 1;
      // Calculate the Pochhammer symbol.
      for(int i = 0; i < j - 1; i++) {
	term *= (0.5 - j + i);
      }
      for(int i = 0; i <= j - 1; i++) {
	sum += term;
	term *= T / (1.5 + i);
      }
      int sign = (j % 2)? 1: -1;
      return (std::tgamma(j + 0.5) * std::erf(std::sqrt(T)) -
	      sign * std::exp(-T) * std::sqrt(T) * sum) /
	(2 * std::pow(T, j + 0.5));
    }
  } else {
    return integrate_square(j, T, this->opts.getintoption("boys points"));
  }
}


static double integrate_linear(int j, double T, int points) {
  double step = 1.0 / (points - 1);
  double sum = 0;
  for(int i = 1; i < points - 1; i++) {
    sum += step * std::pow(i * step, 2 * j) *
      std::exp(-T * (i * step));
  }
  sum += step * std::exp(-T) / 2;
  return sum;
}

static double extrapolant_linear(int i, int k, int j, double T,
			  std::map<std::pair<int, int>, double> &ints) {
  if(ints.count(std::pair<int, int>(i, k)) != 0) {
    return ints.at(std::pair<int, int>(i, k));
  }
  if(i == 0) {
    double res = integrate_linear(j, T, 2 << k);
    ints[std::pair<int, int>(i, k)] = res;
    return res;
  } else {
    double res = extrapolant_linear(i - 1, k + 1, j, T, ints) -
      extrapolant_linear(i - 1, k, j, T, ints) / (std::pow(2, k) - 1);
    ints[std::pair<int, int>(i, k)] = res;
    return res;
  }
}
  

double compchem::AnalyticIntegral::boys_linear(int j, double T) const {
  // If T = 0, this is just a power function.
  if(T == 0) {
    return 1.0 / (2.0 * j + 1.0);
  }
  // If j = 0, this integral is just an exponential function.
  if(j == 0) {
    return (1 - std::exp(-T)) / T;
  }

  if(this->opts.getbooloption("analytic boys")) {
    double sum = 0, coef = 1, pz = 1;
    for(int i = 0; i <= 2 * j; i++) {
      sum += coef * pz;
      coef /= i + 1;
      pz *= T;
    }
    return std::tgamma(2 * j + 1) * (1 - std::exp(-T) * sum) /
      std::pow(T, 2 * j + 1);
  } else {
    return integrate_linear(j, T, this->opts.getintoption("boys points"));
  }
}
