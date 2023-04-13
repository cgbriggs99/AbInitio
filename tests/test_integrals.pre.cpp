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

int test_boys() {
  FILE *fp = fopen("${CMAKE_SOURCE_DIR}/tests/boys_test.dat", "r");

  int warns = 0;
  GlobalOptions::getsingleton().setbooloption("analytic boys", true);
  GlobalOptions::getsingleton().setintoption("boys points", 100);
  GlobalOptions::getsingleton().setintoption("threads", 32);
  
  AnalyticIntegral ints;
  #undef CONV
  #define CONV 1e-5
  while(!feof(fp)) {
    int j;
    double T, res;
    fscanf(fp, "%d\t%lf\t%lf\n", &j, &T, &res);
    double calc = ints.boys_square(j, T);
    ASSERT_WARN_MSG(NEAR(res, calc), warns,
		    "Boys function values are different. Got %lf, expected %lf.\n",
		    calc, res);
  }
  fclose(fp);
  return warns;
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

void read_square(FILE *fp, double *out, int rows) {
  char buf[1024] = {'\0'};
  int swaps[] = {0, 1, 4, 2, 3, 5, 6};

  while(!feof(fp)) {
    fgets(buf, 1023, fp);

    // Read data.
    int mu, nu;
    double val;
    sscanf(buf, "%d %d %lf", &mu, &nu, &val);

    // The files are 1-indexed, but we need zero-indexed.
    mu--;
    nu--;
    // The orbital order is different.
    mu = swaps[mu];
    nu = swaps[nu];
    // File uses symmetry.
    out[mu * rows + nu] = val;
    out[nu * rows + mu] = val;
  }
}

void read_tei(FILE *fp, TEIArray *out) {
  // Clear the output, since any skipped entries are assumed to be zero.
  for(int i = 0; i < out->getsize(); i++) {
    out->at_direct(i) = 0;
  }
  
  char buf[1024] = {'\0'};
  int swaps[] = {0, 1, 4, 2, 3, 5, 6};

  while(!feof(fp)) {
    fgets(buf, 1023, fp);
    
    // Read data.
    int mu, nu, lam, sig;
    double val;
    sscanf(buf, "%d %d %d %d %lf", &mu, &nu, &lam, &sig, &val);

    // File is 1-indexed.
    mu--;
    nu--;
    lam--;
    sig--;

    // Orbitals are in a different order.
    mu = swaps[mu];
    nu = swaps[nu];
    lam = swaps[lam];
    sig = swaps[sig];
    out->at(mu, nu, lam, sig) = val;
  }
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
  
  ASSERT_WARN_MSG(NEAR(boys, exact), warns,
		  "Expected Boys function %lf, got %lf.\n",
		  exact, boys);

  // Generate the integrals.
  factory.Smatrix(mol, S, &sdim);
  factory.Tmatrix(mol, T, &tdim);
  factory.Vmatrix(mol, V, &vdim);
  TEIArray *tei = factory.TEIints(mol);

  // Make sure they at least have the right size.
  ASSERT_WARN_MSG(sdim == 7, warns, "Got %d, expected 7.\n", sdim);
  ASSERT_WARN_MSG(tdim == 7, warns, "Got %d, expected 7.\n", tdim);
  ASSERT_WARN_MSG(tdim == 7, warns, "Got %d, expected 7.\n", tdim);
  ASSERT_WARN_MSG(tei->getdim() == 7, warns, "Got %d, expected 7.\n",
  		  tei->getdim());


  // Read in the pre-computed arrays.
  double *S_pre = new double[49],
    *T_pre = new double[49],
    *V_pre = new double[49];
  TEIArray *tei_pre =  new TEIArray(7);

  FILE *sfile = fopen("${CMAKE_SOURCE_DIR}/tests/S.dat", "r");
  read_square(sfile, S_pre, 7);
  fclose(sfile);

  FILE *tfile = fopen("${CMAKE_SOURCE_DIR}/tests/T.dat", "r");
  read_square(tfile, T_pre, 7);
  fclose(tfile);

  FILE *vfile = fopen("${CMAKE_SOURCE_DIR}/tests/V.dat", "r");
  read_square(vfile, V_pre, 7);
  fclose(vfile);

  FILE *tei_file = fopen("${CMAKE_SOURCE_DIR}/tests/ERI.dat", "r");
  read_tei(tei_file, tei_pre);
  fclose(tei_file);

  // Compare the real vs the calculated.
  for(int i = 0; i < 7; i++) {
    for(int j = 0; j <= i; j++) {
      ASSERT_WARN_MSG(NEAR(S[i * 7 + j], S_pre[i * 7 + j]), warns,
		      "Overlap matrix not the same at (%d, %d). Got %lf, expected %lf.\n",
		      i, j, S[i * 7 + j], S_pre[i * 7 + j]);
      ASSERT_WARN_MSG(NEAR(T[i * 7 + j], T_pre[i * 7 + j]), warns,
		      "Kinetic energy matrix not the same at (%d, %d). Got %lf, expected %lf.\n",
		      i, j, T[i * 7 + j], T_pre[i * 7 + j]);
      ASSERT_WARN_MSG(NEAR(V[i * 7 + j], V_pre[i * 7 + j]), warns,
		      "Potential energy matrix not the same at (%d, %d). Got %lf, expected %lf.\n",
		      i, j, V[i * 7 + j], V_pre[i * 7 + j]);
    }
  }

  for(int i = 0; i < tei->getsize(); i++) {
    int mu, nu, lam, sig;
    tei->indextoquad(i, &mu, &nu, &lam, &sig);
    ASSERT_WARN_MSG(NEAR(tei->at(mu, nu, lam, sig), tei_pre->at(mu, nu, lam, sig)),
		    warns,
		    "Electron repulsion integrals not the same at (%d %d | %d %d). Got %lf, expected %lf.\n",
		    mu, nu, lam, sig, tei->at(mu, nu, lam, sig),
		    tei_pre->at(mu, nu, lam, sig));
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
  delete[] S_pre;
  delete[] T_pre;
  delete[] V_pre;
  delete tei_pre;
  

  return warns;
}

int main(void) {
  int warns = 0, errs = 0, ret;

  ret = test_boys();
  if(ret == -1) {
    errs++;
  } else {
    warns += ret;
  }
  
  ret = test_tei_array();
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
  
  
