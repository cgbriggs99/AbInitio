#include "atom.hpp"
#include <string>

using namespace compchem;

Atom::Atom() {
  this->Z = 0;
  this->charge = 0;
  this->x = 0;
  this->y = 0;
  this->z = 0;
  this->mass = 0;
  this->norbitals = 0;
  this->orbitals = std::vector<BasisOrbital *>(0);
}

Atom::Atom(int Z, double x, double y, double z) {
  this->Z = Z;
  this->charge = 0;
  this->x = x;
  this->y = y;
  this->z = z;
  this->mass = getAbundantMass(Z);
  this->norbitals = 0;
  this->orbitals = std::vector<BasisOrbital *>(0);
}

Atom::Atom(int Z, double mass, double x, double y, double z) {
  this->Z = Z;
  this->charge = 0;
  this->x = x;
  this->y = y;
  this->z = z;
  this->mass = mass;
  this->norbitals = 0;
  this->orbitals = std::vector<BasisOrbital *>(0);
}

Atom::Atom(int Z, double x, double y, double z, int norbs, const BasisOrbital *pass) {
  this->Z = Z;
  this->charge = 0;
  this->x = x;
  this->y = y;
  this->z = z;
  this->mass = getAbundantMass(Z);
  this->norbitals = norbs;
  this->orbitals = std::vector<BasisOrbital *>(norbs);

  for(int i = 0; i < norbs; i++) {
    this->orbitals[i] = &pass[i].copy();
  }
}

Atom::Atom(int Z, double mass, double x, double y, double z, int norbs, const BasisOrbital *pass) {
  this->Z = Z;
  this->charge = 0;
  this->x = x;
  this->y = y;
  this->z = z;
  this->mass = mass;
  this->norbitals = norbs;
  this->orbitals = std::vector<BasisOrbital *>(norbs);

  for(int i = 0; i < norbs; i++) {
    this->orbitals[i] = &pass[i].copy();
  }
}

Atom::Atom(int Z, int charge, double x, double y, double z) {
  this->Z = Z;
  this->charge = charge;
  this->x = x;
  this->y = y;
  this->z = z;
  this->mass = getAbundantMass(Z);
  this->norbitals = 0;
  this->orbitals = std::vector<BasisOrbital *>(0);
}

Atom::Atom(int Z, int charge, double mass, double x, double y, double z) {
  this->Z = Z;
  this->charge = charge;
  this->x = x;
  this->y = y;
  this->z = z;
  this->mass = mass;
  this->norbitals = 0;
  this->orbitals = std::vector<BasisOrbital *>(0);
}

Atom::Atom(int Z, int charge, double x, double y, double z, int norbs, const BasisOrbital *pass) {
  this->Z = Z;
  this->charge = charge;
  this->x = x;
  this->y = y;
  this->z = z;
  this->mass = getAbundantMass(Z);
  this->norbitals = norbs;
  this->orbitals = std::vector<BasisOrbital *>(norbs);

  for(int i = 0; i < norbs; i++) {
    this->orbitals[i] = &pass[i].copy();
  }
}

Atom::Atom(int Z, int charge, double mass, double x, double y, double z, int norbs, const BasisOrbital *pass) {
  this->Z = Z;
  this->charge = charge;
  this->x = x;
  this->y = y;
  this->z = z;
  this->mass = mass;
  this->norbitals = norbs;
  this->orbitals = std::vector<BasisOrbital *>(norbs);

  for(int i = 0; i < norbs; i++) {
    this->orbitals[i] = &pass[i].copy();
  }
}

Atom::Atom(const Atom &copy) {
  this->Z = copy.Z;
  this->charge = copy.charge;
  this->x = copy.x;
  this->y = copy.y;
  this->z = copy.z;
  this->norbitals = copy.norbitals;
  this->orbitals = std::vector<BasisOrbital *>(this->norbitals);

  for(int i = 0; i < this->norbitals; i++) {
    this->orbitals[i] = &(copy.orbitals[i]->copy());
  }
}

Atom::~Atom() {
  for(int i = 0; i < this->norbitals; i++) {
    if(this->orbitals[i] != nullptr) {
      delete this->orbitals[i];
    }
  }
}

double Atom::getx() const {
  return this->x;
}

double Atom::gety() const {
  return this->y;
}

double Atom::getz() const {
  return this->z;
}

double Atom::getmass() const {
  return this->mass;
}

int Atom::getZ() const {
  return this->Z;
}

int Atom::getcharge() const {
  return this->charge;
}

int Atom::getnorbitals() const {
  return this->norbitals;
}

const std::vector<BasisOrbital *> &Atom::getorbitals() const {
  return this->orbitals;
}

int getZFromSymb(const std::string &symb) {
  static std::string symbs[] = {
"H",                                                                                                                                                       "He",
"Li","Be",                                                                                                                        "B", "C", "N", "O", "F", "Ne",
"Na","Mg",                                                                                                                        "Al","Si","P", "S", "Cl","Ar",
"K", "Ca",                                                                      "Sc","Ti","V", "Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
"Rb","Sr",                                                                      "Y", "Zr","Nb","Tc","Mo","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I", "Xe",
"Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W", "Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
"Fr","Ra","Ac","Th","Pa","U", "Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Ru","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"
  };

  // Deuterium and tritium get special symbols.
  if(symb == "D" || symb == "T") {
    return 1;
  }
  for(int i = 0; i < 118; i++) {
    if(symb == symbs[i]) {
      return i + 1;
    }
  }
  return 0;
}

double getAbundantMass(int Z) {
  // Obtained from https://www.ciaaw.org.
  static double masses[] = {
    1.0078250322,
    4.0026032545,
    7.01600344,
    9.0121831,
    11.00930517,
    12, // Exact
    14.003074004,
    15.994914619,
    18.998403162,
    19.99244018,
    22.98976928,
    23.98504170,
    26.9815384,
    27.976926535,
    30.973761998,
    31.972071174,
    34.9688527,
    39.96238312,
    38.96370649,
    39.9625909,
    44.955907,
    47.9479409,
    50.943957,
    51.940505,
    54.938043,
    55.934936,
    58.933194,
    57.935342,
    62.929597,
    63.929142, // Most abundant does not have majority.
    68.925573,
    73.92117776,
    74.921595,
    79.916522,
    78.918338,
    83.91149773,
    84.91178974,
    87.90561226,
    88.905838,
    89.9046988,
    92.90637,
    97.905404, // Extremely even spread, chose Mo98
    97,
    101.904340,
    102.90549,
    105.903480,
    106.90509,
    113.903365,
    114.90387877,
    119.902202,
    120.90381,
    129.90622275,
    126.90447,
    131.90415509,
    132.90545196,
    137.905247,
    138.90636,
    139.90545,
    140.90766,
    141.90773,
    145,
    151.919739,
    152.921237,
    157.924112,
    158.925354,
    163.929181,
    164.930329,
    165.930299,
    168.934219,
    173.93886755,
    174.940777,
    179.94656,
    180.94800,
    183.950933,
    186.955752,
    191.96148,
    192.962924,
    194.964794,
    196.966570,
    201.970644,
    204.974427,
    207.976652,
    208.98040,
    209,
    210,
    222,
    223,
    226,
    227,
    232.03805,
    231.03588,
    238.05079,
    237,
    244,
    243,
    247,
    247,
    251,
    252,
    257,
    258,
    259,
    266,
    267,
    268,
    269,
    270,
    269,
    278,
    281,
    282,
    285,
    286,
    289,
    290,
    293,
    294,
    294};
  if(Z > 0 && Z <= 118) {
    return masses[Z - 1];
  }
  return 0;
}

Atom &Atom::copy() {
  Atom *out = new Atom(*this);
  return *out;
}
