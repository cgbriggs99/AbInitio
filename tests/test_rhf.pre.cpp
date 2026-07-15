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
    *V = new double[49],
    *mux = new double[49],
    *muy = new double[49],
    *muz = new double[49];
  int sdim, tdim, vdim, ddim;

  IntegralFactory<AnalyticIntegral> factory;

  factory.Smatrix(mol, S, &sdim);
  factory.Tmatrix(mol, T, &tdim);
  factory.Vmatrix(mol, V, &vdim);
  factory.dipole(mol, mux, muy, muz, &ddim);
  TEIArray *tei = factory.TEIints(mol);

  ASSERT_FATAL(mol->getelectrons() == 10);

  RHFWfn *wfn = new RHFWfn(mol->getelectrons(), orbs,
			   S, T, V, tei, mux, muy, muz);

  RHF rhf = RHF();

  double energy = rhf.energy(mol, wfn, wfn);

  auto dipole = rhf.dipole(mol, wfn);

  ASSERT_WARN_MSG(NEAR(energy, -74.942079928192), warns,
		  "Expected -74.942079928192, got %lf.\n",
		  energy);

  ASSERT_WARN_MSG(NEAR(dipole[0], 0), warns,
		  "Expected 0, got %lf.\n",
		  dipole[0]);
  ASSERT_WARN_MSG(NEAR(dipole[1], 0), warns,
		  "Expected 0.603521296525, got %lf.\n",
		  dipole[1]);
  ASSERT_WARN_MSG(NEAR(dipole[2], 0), warns,
		  "Expected 0, got %lf.\n",
		  dipole[2]);

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
  
