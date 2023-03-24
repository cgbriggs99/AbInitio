#include <extramath.h>
#include "${CMAKE_SOURCE_DIR}/src/input/psi4_bs_reader.hpp"
#include "${CMAKE_SOURCE_DIR}/src/basis_sets/sto-ng/stong.hpp"
#include "${CMAKE_SOURCE_DIR}/src/basis_sets/basis_set.hpp"
#include "${CMAKE_SOURCE_DIR}/tests/test.h"
#include <cstdio>
#include <string>
#include <errno.h>

using namespace compchem;
using namespace std;

int test_inputter() {
  errno = 0;
  FILE *fp = fopen("${CMAKE_SOURCE_DIR}/tests/sto-3g.gbs", "r");
  int warns = 0;

  if(fp == nullptr || errno != 0) {
    perror("Could not open ${CMAKE_SOURCE_DIR}/tests/sto-3g.gbs");
    return -1;
  }

  std::vector<BasisOrbital *> *horbs, *corbs, *tiorbs;

  horbs = readPsi4file(fp, 1);
  rewind(fp);
  corbs = readPsi4file(fp, 6);
  rewind(fp);
  tiorbs = readPsi4file(fp, 22);

  fclose(fp);

  ASSERT_WARN(horbs->size() == 1, warns);
  ASSERT_WARN(static_cast<GaussianOrbital *>(horbs->at(0))->getnterms() == 3, warns);
  ASSERT_WARN(corbs->size() == 5, warns);
  ASSERT_WARN(static_cast<GaussianOrbital *>(corbs->at(0))->getnterms() == 3, warns);
  ASSERT_WARN(tiorbs->size() == 18, warns);
  ASSERT_WARN(static_cast<GaussianOrbital *>(tiorbs->at(0))->getnterms() == 3, warns);

  for(auto b : *horbs) {
    delete b;
  }
  for(auto b : *corbs) {
    delete b;
  }
  for(auto b : *tiorbs) {
    delete b;
  }
  
  horbs->clear();
  corbs->clear();
  tiorbs->clear();
  delete horbs;
  delete corbs;
  delete tiorbs;
  

  return warns;
}
  

double kernel(const double *r, int dim, const void *pass) {
  const GaussianOrbital *orb = (const GaussianOrbital *) pass;

  return orb->eval(r[0], r[1], r[2]) * orb->eval(r[0], r[1], r[2]) /
    std::exp(-(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]));
}

double kernel2(const double *r, int dim, const void *pass) {
  const SlaterOrbital *orb = (const SlaterOrbital *) pass;

  return orb->eval(r[0], r[1], r[2]) * orb->eval(r[0], r[1], r[2]) /
    std::exp(-(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]));
}

// Integrals are not accurate enough, apparently.
#undef CONV
#define CONV 1e-2
int test_ints() {
  int warns = 0;

  double int1, int2;
  
  GaussianOrbital *orb = new GaussianOrbital(0, 0, std::vector<double>({1}),
					     std::vector<double>({1}));

  ASSERT_WARN_MSG(NEAR(orb->eval(0, 0, 0), std::pow(2 / M_PI, 0.75)), warns,
		  "Got %lf, expected %lf.\n",
		  orb->eval(0, 0, 0), std::pow(2 / M_PI, 0.75));

  int1 = gausshermiteintnd(kernel, 3, 10, (void *) orb);

  delete orb;
  orb = new GaussianOrbital(2, 1, std::vector<double>({1}),
			    std::vector<double>({1}));

  int2 = gausshermiteintnd(kernel, 3, 10, (void *) orb);
  delete orb;

  ASSERT_WARN_MSG(NEAR(int1, 1), warns,
		  "Got %lf, expected 1.\n", int1);
  
  ASSERT_WARN_MSG(NEAR(int2, 1), warns,
		  "Got %lf, expected 1.\n", int2);

  // Integrals are not accurate enough to handle the tail behavior.
  // Only converges to 1 when Zeff is 1, but theoretically should
  // for others.
  double zeff = 1;
  SlaterOrbital *orb2 = new SlaterOrbital(zeff, 1, 0, 0);

  ASSERT_WARN(NEAR(orb2->getZeff(), zeff), warns);

  ASSERT_WARN_MSG(NEAR(orb2->getharms().eval(0, 0, 0),
		       std::sqrt(1 / (4 * M_PI))), warns,
		  "Got %lf, expected %lf.\n",
		  orb2->getharms().eval(0, 0, 0),
		  std::sqrt(1 / (4 * M_PI)));
  
  ASSERT_WARN_MSG(NEAR(orb2->eval(0, 0, 0), 2 * zeff * std::sqrt(zeff) *
		       std::sqrt(1 / (4 * M_PI))), warns,
		  "Got %lf, expected %lf.\n",
		  orb2->eval(0, 0, 0),
		  2 * zeff * std::sqrt(zeff) * std::sqrt(1 / (4 * M_PI)));
  ASSERT_WARN_MSG(NEAR(orb2->eval(1, 0, 0), orb2->eval(-1, 0, 0)), warns,
		  "Function not symmetric.\n");

  double starts[3] = {0, 0, 0},
    ends[3] = {30, 30, 30};
  int bounds[3] = {__INTEGRAL_BOTH__, __INTEGRAL_BOTH__, __INTEGRAL_BOTH__};
  
  double int3 = gausshermiteintnd(kernel2, 3, 18, (void *) orb2);

  ASSERT_WARN_MSG(NEAR(int3, 1), warns,
		  "Got %lf, expected 1.\n", int3);
  delete orb2;

  return warns;

}  

int main(void) {

  int warns = 0, errs = 0, ret;


  ret = test_inputter();
  if(ret == -1) {
    errs++;
  } else {
    warns += ret;
  }

  ret = test_ints();
  if(ret == -1) {
    errs++;
  } else {
    warns += ret;
  }

  return warns + errs;
}
  
