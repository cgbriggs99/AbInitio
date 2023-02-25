#include "../src/util/polynomial.hpp"
#include "test.h"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <errno.h>

using namespace compchem;
using namespace std;

int test_poly() {
  int warns = 0;
  int pow1[] = {1, 0, 0},
    pow2[] = {1, 1, 1},
    pow3[] = {2, 2, 2, 1, 1, 0, 0, 0, 0},
    pow4[] = {1, 1, 0, 2, 2, 2, 0, 0, 0};
  double coef1[] = {1},
    coef2[] = {3},
    coef3[] = {1, -1, 4},
    coef4[] = {-1, 1, 4};
  
  Polynomial<3> *p1 = new Polynomial<3>(10),
    *p2 = new Polynomial<3>(0),
    *p3 = new Polynomial<3>(pow1, coef1, 1),
    *p4 = new Polynomial<3>(pow2, coef2, 1),
    *p5 = new Polynomial<3>(pow3, coef3, 3),
    *p5_copy = new Polynomial<3>(*p5),
    *p6 = new Polynomial<3>(pow4, coef4, 3);

  ASSERT_WARN(*p6 == *p5, warns);
  ASSERT_WARN(*p5 == *p6, warns);
  delete p6;

  // Check to make sure the sizes are right.
  ASSERT_WARN(p1->getsize() == 1, warns);
  ASSERT_WARN(p2->getsize() == 0, warns);
  ASSERT_WARN(p3->getsize() == 1, warns);
  ASSERT_WARN(p4->getsize() == 1, warns);
  ASSERT_WARN(p5->getsize() == 3, warns);
  ASSERT_WARN(p5_copy->getsize() == 3, warns);

  // Check to make sure the coefficients are right.
  ASSERT_WARN(NEAR(p1->getcoef(0), 10), warns);
  ASSERT_WARN(NEAR(p2->getcoef(0), 0), warns);
  ASSERT_WARN(NEAR(p3->getcoef(0), coef1[0]),warns);
  ASSERT_WARN(NEAR(p4->getcoef(0), coef2[0]),warns);
  for(int i = 0; i < 3; i++) {
    ASSERT_WARN(NEAR(p5->getcoef(i), coef3[i]), warns);
    ASSERT_WARN(NEAR(p5_copy->getcoef(i), p5->getcoef(i)), warns);
  }

  // Check to make sure the term orders are right.
  ASSERT_WARN(p2->gettermorder(0) == nullptr, warns);
  for(int i = 0; i < 3; i++) {
    ASSERT_WARN(p1->gettermorder(0)[i] == 0, warns);
    ASSERT_WARN(p3->gettermorder(0)[i] == pow1[i], warns);
    ASSERT_WARN(p4->gettermorder(0)[i] == pow2[i], warns);
    for(int j = 0; j < 3; j++) {
      ASSERT_WARN(p5->gettermorder(j)[i] == pow3[j * 3 + i], warns);
      ASSERT_WARN(p5_copy->gettermorder(j)[i] == p5->gettermorder(j)[i], warns);
    }
  }

  // Check that evaluating works.
  ASSERT_WARN(NEAR(p5->eval(0.0, 0.0, 0.0), 4), warns);
  ASSERT_WARN(NEAR(p2->eval(0.0, 0.0, 0.0), 0), warns);
  ASSERT_WARN(NEAR(p5->eval(2.0, 2.0, 2.0), 64), warns);

  // Check that arithmetic doesn't break.
  int pc1[] = {2, 2, 2, 1, 1, 0, 0, 0, 0};
  double cc1[] = {1, -1, 4};
  Polynomial<3> *check = new Polynomial<3>(pc1, cc1, 3);
  ASSERT_WARN(*p5 == *check, warns);
  if(*p5 != *check) {
    fprintf(stderr, "p5 coefs: [ ");
    for(int i = 0; i < p5->getsize(); i++) {
      fprintf(stderr, "%lf ", p5->getcoef(i));
    }
    fprintf(stderr, "] p5 expos: [ ");
    for(int i = 0; i < p5->getsize() * 3; i++) {
      fprintf(stderr, "%d ", p5->gettermorder(0)[i]);
    }
    fprintf(stderr, "] expected coefs: [ ");
    for(int i = 0; i < check->getsize(); i++) {
      fprintf(stderr, "%lf ", check->getcoef(i));
    }
    fprintf(stderr, "] expected expos: [ ");
    for(int i = 0; i < check->getsize() * 3; i++) {
      fprintf(stderr, "%d ", check->gettermorder(0)[i]);
    }
    perror("]\n");
  }
  delete check;
  
  int pc2[] = {2, 2, 2, 1, 1, 0, 0, 0, 0};
  double cc2[] = {1, -1, 4};
  check = new Polynomial<3>(pc2, cc2, 3);
  *p5 += *p2;
  ASSERT_WARN(*p5 == *check, warns);
  if(*p5 != *check) {
    fprintf(stderr, "p5 coefs: [ ");
    for(int i = 0; i < p5->getsize(); i++) {
      fprintf(stderr, "%lf ", p5->getcoef(i));
    }
    fprintf(stderr, "] p5 expos: [ ");
    for(int i = 0; i < p5->getsize() * 3; i++) {
      fprintf(stderr, "%d ", p5->gettermorder(0)[i]);
    }
    fprintf(stderr, "] expected coefs: [ ");
    for(int i = 0; i < check->getsize(); i++) {
      fprintf(stderr, "%lf ", check->getcoef(i));
    }
    fprintf(stderr, "] expected expos: [ ");
    for(int i = 0; i < check->getsize() * 3; i++) {
      fprintf(stderr, "%d ", check->gettermorder(0)[i]);
    }
    perror("]\n");
  }
  delete check;

  int pc3[] = {2, 2, 2, 1, 1, 0, 0, 0, 0};
  double cc3[] = {1, -1, 14};
  check = new Polynomial<3>(pc3, cc3, 3);
  *p5 += *p1;
  ASSERT_WARN(*p5 == *check, warns);
  if(*p5 != *check) {
    fprintf(stderr, "p5 coefs: [ ");
    for(int i = 0; i < p5->getsize(); i++) {
      fprintf(stderr, "%lf ", p5->getcoef(i));
    }
    fprintf(stderr, "] p5 expos: [ ");
    for(int i = 0; i < p5->getsize() * 3; i++) {
      fprintf(stderr, "%d ", p5->gettermorder(0)[i]);
    }
    fprintf(stderr, "] expected coefs: [ ");
    for(int i = 0; i < check->getsize(); i++) {
      fprintf(stderr, "%lf ", check->getcoef(i));
    }
    fprintf(stderr, "] expected expos: [ ");
    for(int i = 0; i < check->getsize() * 3; i++) {
      fprintf(stderr, "%d ", check->gettermorder(0)[i]);
    }
    perror("]\n");
  }
  delete check;

  int pc4[] = {3, 3, 3, 2, 2, 1, 1, 1, 1};
  double cc4[] = {3, -3, 42};
  check = new Polynomial<3>(pc4, cc4, 3);
  *p5 *= *p4;
  ASSERT_WARN(*p5 == *check, warns);
  if(*p5 != *check) {
    fprintf(stderr, "p5 coefs: [ ");
    for(int i = 0; i < p5->getsize(); i++) {
      fprintf(stderr, "%lf ", p5->getcoef(i));
    }
    fprintf(stderr, "] p5 expos: [ ");
    for(int i = 0; i < p5->getsize() * 3; i++) {
      fprintf(stderr, "%d ", p5->gettermorder(0)[i]);
    }
    fprintf(stderr, "] expected coefs: [ ");
    for(int i = 0; i < check->getsize(); i++) {
      fprintf(stderr, "%lf ", check->getcoef(i));
    }
    fprintf(stderr, "] expected expos: [ ");
    for(int i = 0; i < check->getsize() * 3; i++) {
      fprintf(stderr, "%d ", check->gettermorder(0)[i]);
    }
    perror("]\n");
  }
  delete check;

  int pc5[] = {3, 3, 3, 2, 2, 1, 1, 1, 1};
  double cc5[] = {1, -1, 14};
  check = new Polynomial<3>(pc5, cc5, 3);
  *p5 /= 3;
  ASSERT_WARN(*p5 == *check, warns);
  if(*p5 != *check) {
    fprintf(stderr, "p5 coefs: [ ");
    for(int i = 0; i < p5->getsize(); i++) {
      fprintf(stderr, "%lf ", p5->getcoef(i));
    }
    fprintf(stderr, "] p5 expos: [ ");
    for(int i = 0; i < p5->getsize() * 3; i++) {
      fprintf(stderr, "%d ", p5->gettermorder(0)[i]);
    }
    fprintf(stderr, "] expected coefs: [ ");
    for(int i = 0; i < check->getsize(); i++) {
      fprintf(stderr, "%lf ", check->getcoef(i));
    }
    fprintf(stderr, "] expected expos: [ ");
    for(int i = 0; i < check->getsize() * 3; i++) {
      fprintf(stderr, "%d ", check->gettermorder(0)[i]);
    }
    perror("]\n");
  }
  delete check;

  int pc6[] = {3, 3, 3, 2, 2, 1, 1, 1, 1, 0, 0, 0};
  double cc6[] = {1, -1, 14, 1};
  check = new Polynomial<3>(pc6, cc6, 4);
  *p5 += 1;
  ASSERT_WARN(*p5 == *check, warns);
  if(*p5 != *check) {
    fprintf(stderr, "p5 coefs: [ ");
    for(int i = 0; i < p5->getsize(); i++) {
      fprintf(stderr, "%lf ", p5->getcoef(i));
    }
    fprintf(stderr, "] p5 expos: [ ");
    for(int i = 0; i < p5->getsize() * 3; i++) {
      fprintf(stderr, "%d ", p5->gettermorder(0)[i]);
    }
    fprintf(stderr, "] expected coefs: [ ");
    for(int i = 0; i < check->getsize(); i++) {
      fprintf(stderr, "%lf ", check->getcoef(i));
    }
    fprintf(stderr, "] expected expos: [ ");
    for(int i = 0; i < check->getsize() * 3; i++) {
      fprintf(stderr, "%d ", check->gettermorder(0)[i]);
    }
    perror("]\n");
  }
  delete check;

  int pc7[] = {3, 3, 3, 2, 2, 1, 1, 1, 1, 0, 0, 0};
  double cc7[] = {1, -1, 14, -1};
  check = new Polynomial<3>(pc7, cc7, 4);
  *p5 -= 2;
  ASSERT_WARN(*p5 == *check, warns);
  if(*p5 != *check) {
    fprintf(stderr, "p5 coefs: [ ");
    for(int i = 0; i < p5->getsize(); i++) {
      fprintf(stderr, "%lf ", p5->getcoef(i));
    }
    fprintf(stderr, "] p5 expos: [ ");
    for(int i = 0; i < p5->getsize() * 3; i++) {
      fprintf(stderr, "%d ", p5->gettermorder(0)[i]);
    }
    fprintf(stderr, "] expected coefs: [ ");
    for(int i = 0; i < check->getsize(); i++) {
      fprintf(stderr, "%lf ", check->getcoef(i));
    }
    fprintf(stderr, "] expected expos: [ ");
    for(int i = 0; i < check->getsize() * 3; i++) {
      fprintf(stderr, "%d ", check->gettermorder(0)[i]);
    }
    perror("]\n");
  }
  delete check;

  int pc8[] = {3, 3, 3, 2, 2, 1, 1, 1, 1, 0, 0, 0};
  double cc8[] = {10, -10, 140, -10};
  check = new Polynomial<3>(pc8, cc8, 4);
  *p5 *= 10;
  ASSERT_WARN(*p5 == *check, warns);
  if(*p5 != *check) {
    fprintf(stderr, "p5 coefs: [ ");
    for(int i = 0; i < p5->getsize(); i++) {
      fprintf(stderr, "%lf ", p5->getcoef(i));
    }
    fprintf(stderr, "] p5 expos: [ ");
    for(int i = 0; i < p5->getsize() * 3; i++) {
      fprintf(stderr, "%d ", p5->gettermorder(0)[i]);
    }
    fprintf(stderr, "] expected coefs: [ ");
    for(int i = 0; i < check->getsize(); i++) {
      fprintf(stderr, "%lf ", check->getcoef(i));
    }
    fprintf(stderr, "] expected expos: [ ");
    for(int i = 0; i < check->getsize() * 3; i++) {
      fprintf(stderr, "%d ", check->gettermorder(0)[i]);
    }
    perror("]\n");
  }
  delete check;

  int pc9[] = {3, 3, 3, 2, 2, 1, 1, 1, 1, 0, 0, 0};
  double cc9[] = {-10, 10, -140, 10};
  check = new Polynomial<3>(pc9, cc9, 4);
  *p5 = -*p5;
  ASSERT_WARN(*p5 == *check, warns);
  if(*p5 != *check) {
    fprintf(stderr, "p5 coefs: [ ");
    for(int i = 0; i < p5->getsize(); i++) {
      fprintf(stderr, "%lf ", p5->getcoef(i));
    }
    fprintf(stderr, "] p5 expos: [ ");
    for(int i = 0; i < p5->getsize() * 3; i++) {
      fprintf(stderr, "%d ", p5->gettermorder(0)[i]);
    }
    fprintf(stderr, "] expected coefs: [ ");
    for(int i = 0; i < check->getsize(); i++) {
      fprintf(stderr, "%lf ", check->getcoef(i));
    }
    fprintf(stderr, "] expected expos: [ ");
    for(int i = 0; i < check->getsize() * 3; i++) {
      fprintf(stderr, "%d ", check->gettermorder(0)[i]);
    }
    perror("]\n");
  }
  delete check;

  int pc10[] = {0, 0, 0};
  double cc10[] = {10};
  check = new Polynomial<3>(pc10, cc10, 1);
  *p1 = 10;
  ASSERT_WARN(*p1 == *check, warns);
  if(*p1 != *check) {
    fprintf(stderr, "p1 coefs: [ ");
    for(int i = 0; i < p1->getsize(); i++) {
      fprintf(stderr, "%lf ", p1->getcoef(i));
    }
    fprintf(stderr, "] p1 expos: [ ");
    for(int i = 0; i < p1->getsize() * 3; i++) {
      fprintf(stderr, "%d ", p1->gettermorder(0)[i]);
    }
    fprintf(stderr, "] expected coefs: [ ");
    for(int i = 0; i < check->getsize(); i++) {
      fprintf(stderr, "%lf ", check->getcoef(i));
    }
    fprintf(stderr, "] expected expos: [ ");
    for(int i = 0; i < check->getsize() * 3; i++) {
      fprintf(stderr, "%d ", check->gettermorder(0)[i]);
    }
    perror("]\n");
  }
  
  *p2 = 10;
  ASSERT_WARN(*p2 == *check, warns);
  if(*p2 != *check) {
    fprintf(stderr, "p2 coefs: [ ");
    for(int i = 0; i < p2->getsize(); i++) {
      fprintf(stderr, "%lf ", p2->getcoef(i));
    }
    fprintf(stderr, "] p2 expos: [ ");
    for(int i = 0; i < p2->getsize() * 3; i++) {
      fprintf(stderr, "%d ", p2->gettermorder(0)[i]);
    }
    fprintf(stderr, "] expected coefs: [ ");
    for(int i = 0; i < check->getsize(); i++) {
      fprintf(stderr, "%lf ", check->getcoef(i));
    }
    fprintf(stderr, "] expected expos: [ ");
    for(int i = 0; i < check->getsize() * 3; i++) {
      fprintf(stderr, "%d ", check->gettermorder(0)[i]);
    }
    perror("]\n");
  }
  delete check;

  *p2 = 0;
  ASSERT_WARN(*p2 == 0, warns);
  if(*p2 != 0) {
    fprintf(stderr, "p2 coefs: [ ");
    for(int i = 0; i < p2->getsize(); i++) {
      fprintf(stderr, "%lf ", p2->getcoef(i));
    }
    fprintf(stderr, "] p2 expos: [ ");
    for(int i = 0; i < p2->getsize() * 3; i++) {
      fprintf(stderr, "%d ", p2->gettermorder(0)[i]);
    }
    perror("]\n");
  }

  // Check that comparisons don't break.
  ASSERT_WARN(*p1 == 10, warns);
  ASSERT_WARN(*p2 == 0, warns);
  ASSERT_WARN(!(*p1 == *p2), warns);
  ASSERT_WARN(*p1 == *p1, warns);
  ASSERT_WARN(*p2 == 0, warns);
  ASSERT_WARN(*p5 != *p5_copy, warns);
  ASSERT_WARN(!(*p3 != *p3), warns);

  // Check that more arithmetic doesn't break.
  int pc11[] = {2, 0, 0, 0, 0, 0};
  double cc11[] = {2, 10};
  check = new Polynomial<3>(pc11, cc11, 2);
  *p1 += 2 * compchem::pow(*p3, 2);
  ASSERT_WARN(*p1 == *check, warns);
  if(*p1 != *check) {
    fprintf(stderr, "p1 coefs: [ ");
    for(int i = 0; i < p1->getsize(); i++) {
      fprintf(stderr, "%lf ", p1->getcoef(i));
    }
    fprintf(stderr, "] p1 expos: [ ");
    for(int i = 0; i < p1->getsize() * 3; i++) {
      fprintf(stderr, "%d ", p1->gettermorder(0)[i]);
    }
    fprintf(stderr, "] expected coefs: [ ");
    for(int i = 0; i < check->getsize(); i++) {
      fprintf(stderr, "%lf ", check->getcoef(i));
    }
    fprintf(stderr, "] expected expos: [ ");
    for(int i = 0; i < check->getsize() * 3; i++) {
      fprintf(stderr, "%d ", check->gettermorder(0)[i]);
    }
    perror("]\n");
  }
  delete check;

  int pc12[] = {2, 0, 0, 0, 0, 0};
  double cc12[] = {2, 10};
  check = new Polynomial<3>(pc12, cc12, 2);
  *p3 = *p1 + *p2;
  ASSERT_WARN(*p3 == *check, warns);
  if(*p3 != *check) {
    fprintf(stderr, "p3 coefs: [ ");
    for(int i = 0; i < p3->getsize(); i++) {
      fprintf(stderr, "%lf ", p3->getcoef(i));
    }
    fprintf(stderr, "] p3 expos: [ ");
    for(int i = 0; i < p3->getsize() * 3; i++) {
      fprintf(stderr, "%d ", p3->gettermorder(0)[i]);
    }
    fprintf(stderr, "] expected coefs: [ ");
    for(int i = 0; i < check->getsize(); i++) {
      fprintf(stderr, "%lf ", check->getcoef(i));
    }
    fprintf(stderr, "] expected expos: [ ");
    for(int i = 0; i < check->getsize() * 3; i++) {
      fprintf(stderr, "%d ", check->gettermorder(0)[i]);
    }
    perror("]\n");
  }
  delete check;

  *p1 = *p2 + 10;
  ASSERT_WARN(*p1 == 10, warns);
  
  *p1 = 10 + *p2;
  ASSERT_WARN(*p1 == 10, warns);

  int pc13[] = {2, 0, 0, 3, 3, 3, 2, 2, 1, 1, 1, 1};
  double cc13[] = {2, 10, -10, 140};
  check = new Polynomial<3>(pc13, cc13, 4);
  *p4 = *p3 - *p5;
  ASSERT_WARN(*p4 == *check, warns);
  if(*p4 != *check) {
    fprintf(stderr, "p4 coefs: [ ");
    for(int i = 0; i < p4->getsize(); i++) {
      fprintf(stderr, "%lf ", p4->getcoef(i));
    }
    fprintf(stderr, "] p4 expos: [ ");
    for(int i = 0; i < p4->getsize() * 3; i++) {
      fprintf(stderr, "%d ", p4->gettermorder(0)[i]);
    }
    fprintf(stderr, "] expected coefs: [ ");
    for(int i = 0; i < check->getsize(); i++) {
      fprintf(stderr, "%lf ", check->getcoef(i));
    }
    fprintf(stderr, "] expected expos: [ ");
    for(int i = 0; i < check->getsize() * 3; i++) {
      fprintf(stderr, "%d ", check->gettermorder(0)[i]);
    }
    perror("]\n");
  }
  delete check;

  int pc14[] = {2, 0, 0, 0, 0, 0};
  double cc14[] = {2, 9};
  check = new Polynomial<3>(pc14, cc14, 2);
  *p4 = *p3 - 1;
  ASSERT_WARN(*p4 == *check, warns);
  if(*p4 != *check) {
    fprintf(stderr, "p4 coefs: [ ");
    for(int i = 0; i < p4->getsize(); i++) {
      fprintf(stderr, "%lf ", p4->getcoef(i));
    }
    fprintf(stderr, "] p4 expos: [ ");
    for(int i = 0; i < p4->getsize() * 3; i++) {
      fprintf(stderr, "%d ", p4->gettermorder(0)[i]);
    }
    fprintf(stderr, "] expected coefs: [ ");
    for(int i = 0; i < check->getsize(); i++) {
      fprintf(stderr, "%lf ", check->getcoef(i));
    }
    fprintf(stderr, "] expected expos: [ ");
    for(int i = 0; i < check->getsize() * 3; i++) {
      fprintf(stderr, "%d ", check->gettermorder(0)[i]);
    }
    perror("]\n");
  }
  delete check;

  int pc15[] = {2, 0, 0, 0, 0, 0};
  double cc15[] = {-2, -9};
  check = new Polynomial<3>(pc15, cc15, 2);
  *p4 = 1 - *p3;
  ASSERT_WARN(*p4 == *check, warns);
  if(*p4 != *check) {
    fprintf(stderr, "p4 coefs: [ ");
    for(int i = 0; i < p4->getsize(); i++) {
      fprintf(stderr, "%lf ", p4->getcoef(i));
    }
    fprintf(stderr, "] p4 expos: [ ");
    for(int i = 0; i < p4->getsize() * 3; i++) {
      fprintf(stderr, "%d ", p4->gettermorder(0)[i]);
    }
    fprintf(stderr, "] expected coefs: [ ");
    for(int i = 0; i < check->getsize(); i++) {
      fprintf(stderr, "%lf ", check->getcoef(i));
    }
    fprintf(stderr, "] expected expos: [ ");
    for(int i = 0; i < check->getsize() * 3; i++) {
      fprintf(stderr, "%d ", check->gettermorder(0)[i]);
    }
    perror("]\n");
  }
  delete check;

  int pc16[] = {2, 0, 0, 0, 0, 0};
  double cc16[] = {-20, -90};
  check = new Polynomial<3>(pc16, cc16, 2);
  *p3 = *p1 * *p4;
  ASSERT_WARN(*p3 == *check, warns);
  if(*p3 != *check) {
    fprintf(stderr, "p3 coefs: [ ");
    for(int i = 0; i < p3->getsize(); i++) {
      fprintf(stderr, "%lf ", p3->getcoef(i));
    }
    fprintf(stderr, "] p3 expos: [ ");
    for(int i = 0; i < p3->getsize() * 3; i++) {
      fprintf(stderr, "%d ", p3->gettermorder(0)[i]);
    }
    fprintf(stderr, "] expected coefs: [ ");
    for(int i = 0; i < check->getsize(); i++) {
      fprintf(stderr, "%lf ", check->getcoef(i));
    }
    fprintf(stderr, "] expected expos: [ ");
    for(int i = 0; i < check->getsize() * 3; i++) {
      fprintf(stderr, "%d ", check->gettermorder(0)[i]);
    }
    perror("]\n");
  }
  delete check;
  
  *p4 = *p3 * *p2;
  ASSERT_WARN(*p4 == 0, warns);

  int pc17[] = {2, 0, 0, 0, 0, 0};
  double cc17[] = {-200, -900};
  check = new Polynomial<3>(pc17, cc17, 2);
  *p4 = *p3 * 10;
  ASSERT_WARN(*p4 == *check, warns);
  if(*p4 != *check) {
    fprintf(stderr, "p4 coefs: [ ");
    for(int i = 0; i < p4->getsize(); i++) {
      fprintf(stderr, "%lf ", p4->getcoef(i));
    }
    fprintf(stderr, "] p4 expos: [ ");
    for(int i = 0; i < p4->getsize() * 3; i++) {
      fprintf(stderr, "%d ", p4->gettermorder(0)[i]);
    }
    fprintf(stderr, "] expected coefs: [ ");
    for(int i = 0; i < check->getsize(); i++) {
      fprintf(stderr, "%lf ", check->getcoef(i));
    }
    fprintf(stderr, "] expected expos: [ ");
    for(int i = 0; i < check->getsize() * 3; i++) {
      fprintf(stderr, "%d ", check->gettermorder(0)[i]);
    }
    perror("]\n");
  }
  delete check;

  int pc18[] = {0, 0, 0};
  double cc18[] = {20};
  check = new Polynomial<3>(pc18, cc18, 1);
  *p3 = *p1 * 2;
  ASSERT_WARN(*p3 == *check, warns);
  if(*p3 != *check) {
    fprintf(stderr, "p3 coefs: [ ");
    for(int i = 0; i < p3->getsize(); i++) {
      fprintf(stderr, "%lf ", p3->getcoef(i));
    }
    fprintf(stderr, "] p3 expos: [ ");
    for(int i = 0; i < p3->getsize() * 3; i++) {
      fprintf(stderr, "%d ", p3->gettermorder(0)[i]);
    }
    fprintf(stderr, "] expected coefs: [ ");
    for(int i = 0; i < check->getsize(); i++) {
      fprintf(stderr, "%lf ", check->getcoef(i));
    }
    fprintf(stderr, "] expected expos: [ ");
    for(int i = 0; i < check->getsize() * 3; i++) {
      fprintf(stderr, "%d ", check->gettermorder(0)[i]);
    }
    perror("]\n");
  }
  delete check;
  
  *p2 = *p2 * 3;
  ASSERT_WARN(*p2 == 0, warns);

  int pc19[] = {2, 0, 0, 0, 0, 0};
  double cc19[] = {200.0/3, 310};
  check = new Polynomial<3>(pc19, cc19, 2);
  *p5 = *p1 + *p2 * *p3 - *p4 / 3;
  ASSERT_WARN(*p5 == *check, warns);
  if(*p5 != *check) {
    fprintf(stderr, "p5 coefs: [ ");
    for(int i = 0; i < p5->getsize(); i++) {
      fprintf(stderr, "%lf ", p5->getcoef(i));
    }
    fprintf(stderr, "] p5 expos: [ ");
    for(int i = 0; i < p5->getsize() * 3; i++) {
      fprintf(stderr, "%d ", p5->gettermorder(0)[i]);
    }
    fprintf(stderr, "] expected coefs: [ ");
    for(int i = 0; i < check->getsize(); i++) {
      fprintf(stderr, "%lf ", check->getcoef(i));
    }
    fprintf(stderr, "] expected expos: [ ");
    for(int i = 0; i < check->getsize() * 3; i++) {
      fprintf(stderr, "%d ", check->gettermorder(0)[i]);
    }
    perror("]\n");
  }
  delete check;

  int pc20[] = {4, 0, 0, 2, 0, 0, 0, 0, 0};
  double cc20[] = {80000.0/9, 248000.0/3, 192200};
  check = new Polynomial<3>(pc20, cc20, 3);
  *p1 = 2 * compchem::pow(*p5, 2);
  ASSERT_WARN(p1->almost(*check), warns);
  if(!p1->almost(*check)) {
    fprintf(stderr, "p1 coefs: [ ");
    for(int i = 0; i < p1->getsize(); i++) {
      fprintf(stderr, "%lf ", p1->getcoef(i));
    }
    fprintf(stderr, "] p1 expos: [ ");
    for(int i = 0; i < p1->getsize() * 3; i++) {
      fprintf(stderr, "%d ", p1->gettermorder(0)[i]);
    }
    fprintf(stderr, "] expected coefs: [ ");
    for(int i = 0; i < check->getsize(); i++) {
      fprintf(stderr, "%lf ", check->getcoef(i));
    }
    fprintf(stderr, "] expected expos: [ ");
    for(int i = 0; i < check->getsize() * 3; i++) {
      fprintf(stderr, "%d ", check->gettermorder(0)[i]);
    }
    perror("]\n");
  }
  delete check;
  
  delete p1;
  delete p2;
  delete p3;
  delete p4;
  delete p5;
  delete p5_copy;

  // Test that the spherical harmonics work.
  int pc21[] = {1, 0, 0};
  double cc21[] = {sqrt(3 / (4 * M_PI))};
  check = new Polynomial<3>(pc21, cc21, 1);
  p1 = sphereharm(1, 1);
  ASSERT_WARN(*p1 == *check, warns);
  if(*p1 != *check) {
    fprintf(stderr, "p1 coefs: [ ");
    for(int i = 0; i < p5->getsize(); i++) {
      fprintf(stderr, "%lf ", p5->getcoef(i));
    }
    fprintf(stderr, "] p1 expos: [ ");
    for(int i = 0; i < p1->getsize() * 3; i++) {
      fprintf(stderr, "%d ", p1->gettermorder(0)[i]);
    }
    fprintf(stderr, "] expected coefs: [ ");
    for(int i = 0; i < check->getsize(); i++) {
      fprintf(stderr, "%lf ", check->getcoef(i));
    }
    fprintf(stderr, "] expected expos: [ ");
    for(int i = 0; i < check->getsize() * 3; i++) {
      fprintf(stderr, "%d ", check->gettermorder(0)[i]);
    }
    perror("]\n");
  }
  delete check;
  delete p1;
  
  return warns;
}

int main(void) {
  int warns, errs = 0;

  int ret = test_poly();
  if(ret == -1) {
    errs++;
  } else {
    warns += ret;
  }

  return warns + errs;
}
  

  
