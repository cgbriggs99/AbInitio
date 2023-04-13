#include "${CMAKE_SOURCE_DIR}/src/scf/rhf.hpp"
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

using namespace compchem;

int test_rhf(void) {
  int warns = 0, orbs;

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
  int sdim, tdim, vdim;

  IntegralFactory<AnalyticIntegral> factory;

  factory.Smatrix(mol, S, &sdim);
  factory.Tmatrix(mol, T, &tdim);
  factory.Vmatrix(mol, V, &vdim);
  TEIArray *tei = factory.TEIints(mol);

  ASSERT_FATAL(mol->getelectrons() == 10);

  RHFWfn *wfn = new RHFWfn(mol->getelectrons(), orbs,
			   S, T, V, tei);

  RHF rhf = RHF();

  double energy = rhf.energy(mol, wfn, wfn);

  ASSERT_WARN_MSG(NEAR(energy, -74.942079928192), warns,
		  "Expected -74.942079928192, got %lf.\n",
		  energy);

  delete wfn;
  delete mol;

  return warns;
}

int main(void) {
  int warns, errs;

  int ret = test_rhf();
  if(ret == -1) {
    errs++;
  } else {
    warns += ret;
  }

  return warns;
}
  
