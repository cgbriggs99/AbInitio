#include "${CMAKE_SOURCE_DIR}/src/geometric/rotor.hpp"
#include "${CMAKE_SOURCE_DIR}/src/util/molecule.hpp"
#include <stdlib.h>
#include <string.h>
#include <math.h>
// To be replaced with the appropriate lapacke header.
#include <${LAPACKE_HEADER}>
#include <array>
#include "${CMAKE_SOURCE_DIR}/src/util/constants.hpp"


using namespace compchem;

#define MUCH_LESS 100

rotor_t compchem::comprotor(const Molecule &mol, std::array<double, 3> *out) {

  double *moment = new double[9];
  double comx = mol.getcomx(), comy = mol.getcomy(), comz = mol.getcomz();
  double *real = new double[3], *imag = new double[3];

  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      moment[i + 3 * j] = 0;
    }
  }

  for(int k = 0; k < mol.getsize(); k++) {
    double mass = mol.getatom(k).getmass(),
      x = mol.getatom(k).getx(),
      y = mol.getatom(k).gety(),
      z = mol.getatom(k).getz();
    moment[0] += mass * ((y - comy) * (y - comy) + (z - comz) * (z - comz));
    moment[1] += mass * (x - comx) * (y - comy);
    moment[2] += mass * (x - comx) * (z - comz);
    moment[3] += mass * (y - comy) * (x - comx);
    moment[4] += mass * ((x - comx) * (x - comx) + (z - comz) * (z - comz));
    moment[5] += mass * (y - comy) * (z - comz);
    moment[6] += mass * (z - comz) * (x - comx);
    moment[7] += mass * (z - comz) * (y - comy);
    moment[8] += mass * ((x - comx) * (x - comx) + (y - comy) * (y - comy));
  }

  LAPACKE_dgeev(LAPACK_COL_MAJOR, 'N', 'N', 3, moment, 3, real, imag, nullptr, 1, nullptr, 1);

  // Sort the eigenvalues.
  if(real[0] > real[1]) {
    double swap = real[0];
    real[0] = real[1];
    real[1] = swap;
  }
  if(real[1] > real[2]){
    double swap = real[1];
    real[1] = real[2];
    real[2] = swap;
  }
  if(real[0] > real[1]) {
    double swap = real[0];
    real[0] = real[1];
    real[1] = swap;
  }
  // Convert to joules.
  (*out)[0] = PLANCK_SI / (8 * M_PI * M_PI * LIGHT_SI * AU_TO_KG *
			 BOHR_TO_M * BOHR_TO_M * real[0]);
  (*out)[1] = PLANCK_SI / (8 * M_PI * M_PI * LIGHT_SI * AU_TO_KG *
			 BOHR_TO_M * BOHR_TO_M * real[1]);
  (*out)[2] = PLANCK_SI / (8 * M_PI * M_PI * LIGHT_SI * AU_TO_KG *
			 BOHR_TO_M * BOHR_TO_M * real[2]);

  delete[] imag;
  delete[] moment;

  // Determine the rotor type.
  if(real[0] == real[1] && real[1] == real[2]) {
    delete[] real;
    return SPHERICAL;
  } else if(real[1] == real[2] && real[0] < real[1] / MUCH_LESS) {
    delete[] real;
    return LINEAR;
  } else if(real[0] == real[1] && real[1] != real[2]) {
    delete[] real;
    return OBLATE;
  } else if(real[1] == real[2] && real[0] != real[1]) {
    delete[] real;
    return PROLATE;
  } else {
    delete[] real;
    return ASSYMETRIC;
  }
}
