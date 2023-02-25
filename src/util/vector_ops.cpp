#include "vector_ops.hpp"
#include <vector>
#include <cmath>

using namespace compchem;
using namespace std;

double compchem::rmsnorm(const std::vector<double> &v1,
			 const std::vector<double> &v2) {
  if(v1.size() != v2.size()) {
    return NAN;
  }
  double sum = 0;

  for(int i = 0; i < v1.size(); i++) {
    sum += (v1[i] - v2[i]) * (v1[i] - v2[i]);
  }
  return sqrt(sum / v1.size());
}

double compchem::dotprod(const std::vector<double> &v1,
			 const std::vector<double> &v2) {
  if(v1.size() != v2.size()) {
    return NAN;
  }

  double sum = 0;
  for(int i = 0; i < v1.size(); i++) {
    sum += v1[i] * v2[i];
  }
  return sum;
}

double compchem::vecmag(const std::vector<double> &v1) {
  double sum = 0;
  for(int i = 0; i < v1.size(); i++) {
    sum += v1[i] * v1[i];
  }
  return sqrt(sum);
}
    
