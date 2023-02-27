#ifndef __VECTOR_OPS_HPP__
#define __VECTOR_OPS_HPP__

#include <vector>
namespace compchem {
  
  // Norms.
  double rmsnorm(const std::vector<double> &v1, const std::vector<double> &v2);

  // Dot product
  double dotprod(const std::vector<double> &v1, const std::vector<double> &v2);

  // Magnitude
  double vecmag(const std::vector<double> &vec);

}


#endif
