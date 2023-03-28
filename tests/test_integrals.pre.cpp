#include <extramath.h>
#include "${CMAKE_SOURCE_DIR}/src/input/psi4_bs_reader.hpp"
#include "${CMAKE_SOURCE_DIR}/src/basis_sets/sto-ng/stong.hpp"
#include "${CMAKE_SOURCE_DIR}/src/basis_sets/basis_set.hpp"
#include "${CMAKE_SOURCE_DIR}/tests/test.h"
#include "${CMAKE_SOURCE_DIR}/src/integrals/integrals.hpp"
#include "${CMAKE_SOURCE_DIR}/src/input/input.hpp"
#include "${CMAKE_SOURCE_DIR}/src/opts/default_options.hpp"
#include <cstdio>
#include <string>
#include <errno.h>

using namespace std;
using namespace compchem;

void print_matrix(const double *mat, int rc, const string &name,
		  FILE *fp = stdout) {
  fprintf(fp, name.c_str());
  fprintf(fp, ": \n");
  for(int i = 0; i < rc; i++) {
    fprintf(fp, "[ ");
    for(int j = 0; j < rc; j++) {
      fprintf(fp, "%0.5lf ", mat[i * rc + j]);
    }
    fprintf(fp, "]\n");
  }
}

void print_tei(const TEIArray &tei, const string &name,
	       FILE *fp = stdout) {
  fprintf(fp, name.c_str());
  fprintf(fp, ": \n");
  for(int i = 0; i < tei.getsize(); i++) {
    int mu, nu, lam, sig;
    tei.indextoquad(i, &mu, &nu, &lam, &sig);
    fprintf(fp, "%d %d %d %d: %0.5lf\n", mu, nu, lam, sig, tei.getdata()[i]);
  }
}

int test_tei_array() {
  TEIArray *tei = new TEIArray(7);

  int mu, nu, lam, sig;
  int warns = 0;

  for(int i = 0; i < tei->getsize(); i++) {
    tei->indextoquad(i, &mu, &nu, &lam, &sig);
    ASSERT_WARN_MSG(mu >= 0 && mu < 7 &&
		    nu >= 0 && nu < 7 &&
		    lam >= 0 && lam < 7 &&
		    sig >= 0 && sig < 7, warns,
		    "Index did not work! %d gives %d, %d, %d, %d.\n",
		    i, mu, nu, lam, sig);
  }

  delete tei;

  return warns;
}

int test_integrals() {
  int warns = 0, orbs = 0;
  FILE *geo = fopen("${CMAKE_SOURCE_DIR}/tests/water2.xyz", "r");
  // generate the molecule.
  Molecule *mol = &parseXYZ(geo);
  mol->translate_to_com();
  fclose(geo);
  FILE *basis = fopen("${CMAKE_SOURCE_DIR}/tests/sto-3g.gbs", "r");
  // Get the basis set.
  for(Atom *a : mol->getatoms()) {
    std::vector<BasisOrbital *> *orb = readPsi4file(basis, a->getZ());
    a->setorbitals(*orb);
    for(auto b : *orb) {
      delete b;
    }
    orb->clear();
    delete orb;
    orbs += a->getnorbitals();
  }
  fclose(basis);
  ASSERT_WARN_MSG(orbs == 7, warns, "Got %d, expected 7.\n", orbs);
  
  double *S = new double[49],
    *T = new double[49],
    *V = new double[49];
  int sdim = 0, tdim = 0, vdim = 0;

  GlobalOptions::getsingleton().setbooloption("analytic boys", true);
  GlobalOptions::getsingleton().setintoption("boys points", 100);
  GlobalOptions::getsingleton().setintoption("threads", 32);
  
  IntegralFactory<AnalyticIntegral> factory;

  // Check that the Boys function is correct.
  double x = 5;
  AnalyticIntegral ints;
  double boys = ints.boys_square(2, x);
  double exact = (3 * std::sqrt(x * M_PI) * std::erf(std::sqrt(x)) -
		  130 / std::exp(5)) / 1000 ,
    gam = exact * 2 * std::pow(x, 1.5);
  //ASSERT_WARN_MSG(NEAR(ingamma(1.5, x), gam), warns,
  //		  "Expected incomplete gamma %lf, got %lf.\n",
  //		  gam, ingamma(1.5, x));
  ASSERT_WARN_MSG(NEAR(boys, exact), warns,
		  "Expected Boys function %lf, got %lf.\n",
		  exact, boys);
  
  factory.Smatrix(mol, S, &sdim);
  factory.Tmatrix(mol, T, &tdim);
  factory.Vmatrix(mol, V, &vdim);
  TEIArray *tei = factory.TEIints(mol);

  ASSERT_WARN_MSG(sdim == 7, warns, "Got %d, expected 7.\n", sdim);
  ASSERT_WARN_MSG(tdim == 7, warns, "Got %d, expected 7.\n", tdim);
  ASSERT_WARN_MSG(tdim == 7, warns, "Got %d, expected 7.\n", tdim);
  ASSERT_WARN_MSG(tei->getdim() == 7, warns, "Got %d, expected 7.\n",
  		  tei->getdim());

  for(int i = 0; i < 7; i++) {
    ASSERT_WARN_MSG(NEAR(S[i * 7 + i], 1), warns,
		    "Diagonal entries are not normailzed: %d: %lf\n",
		    i, S[i * 7 + i]);
  }

  for(int i = 0; i < 49; i++) {
    ASSERT_WARN_MSG(std::fabs(S[i]) <= 1 + CONV, warns,
		    "Not correct at %d, %d: %lf\n", i / 7, i % 7, S[i]);
    ASSERT_WARN_MSG(isfinite(T[i]), warns,
		    "Not finite at %d, %d\n", i / 7, i % 7);
    ASSERT_WARN_MSG(isfinite(V[i]), warns,
		    "Not finite at %d, %d\n", i / 7, i % 7);
  }

  for(int i = 0; i < tei->getsize(); i++) {
    ASSERT_WARN(isfinite(tei->getdata()[i]), warns);
  }

  double *H = new double[49];

  for(int i = 0; i < 49; i++) {
    H[i] = T[i] + V[i];
  }
  
  FILE *mats = fopen("mats.dat", "w+");
  print_matrix(S, 7, "Overlap", mats);
  print_matrix(T, 7, "Kinetic Energy", mats);
  print_matrix(V, 7, "Coulomb", mats);
  print_matrix(H, 7, "Hamiltonian", mats);
  print_tei(*tei, "Two-electron integrals", mats);
  fclose(mats);

  delete[] S;
  delete[] T;
  delete[] V;
  delete[] H;
  delete tei;
  delete mol;

  return warns;
}

int main(void) {
  int warns = 0, errs = 0;

  int ret = test_tei_array();
  if(ret == -1) {
    errs++;
  } else {
    warns += ret;
  }
  
  ret = test_integrals();
  if(ret == -1) {
    errs++;
  } else {
    warns += ret;
  }

  return warns + errs;
}
  
  
