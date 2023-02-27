/*
 * electron_configs.cpp
 *
 *  Created on: Feb 7, 2023
 *      Author: connor
 */

#include "electron_configs.hpp"
#include <cstdlib>
#include <stdexcept>
#include <vector>

using namespace compchem;

GSConfig::GSConfig(int z, std::initializer_list<int> conf) {
	this->Z = z;
	this->size_confs = conf.size();
	this->confs = new int[this->size_confs];
	int i = 0;
	for(int gs : conf) {
		this->confs[i] = gs;
		i++;
	}
}

GSConfig::GSConfig(int z, const std::vector<int> &conf) {
	this->Z = z;
	this->size_confs = conf.size();
	this->confs = new int[this->size_confs];
        for(int i = 0; i < this->size_confs; i++) {
	  this->confs[i] = conf[i];
	}
}

GSConfig::~GSConfig() {
  delete[] this->confs;
}

int GSConfig::getShells() const {
  return this->size_confs;
}

int GSConfig::getZ() const {
	return this->Z;
}

int GSConfig::operator[](int index) const {
	return this->confs[index];
}

GSConfig *compchem::getconfig(int Z) {
  // Try to find a way to do this without allocating every time.
  std::vector<int> gsconfigs[] = {
    std::vector<int>({1}),	// H: 1s1
    std::vector<int>({2}),	// He: 1s2
    std::vector<int>({2, 1}),	// Li: 1s2 2s1
    std::vector<int>({2, 2}),	// Be: 1s2 2s2
    std::vector<int>({2, 2, 1}),	// B: 1s2 2s2 2p1
    std::vector<int>({2, 2, 2}),	// C: 1s2 2s2 2p2
    std::vector<int>({2, 2, 3}),	// N: 1s2 2s2 2p3
    std::vector<int>({2, 2, 4}),	// O: 1s2 2s2 2p4
    std::vector<int>({2, 2, 5}),	// F: 1s2 2s2 2p5
    std::vector<int>({2, 2, 6}),	// Ne: 1s2 2s2 2p6
    std::vector<int>({2, 2, 6, 1}),	// Na: [Ne] 3s1
    std::vector<int>({2, 2, 6, 2}),	// Mg: [Ne] 3s2
    std::vector<int>({2, 2, 6, 2, 1}),	// Al: [Ne] 3s2 3p1
    std::vector<int>({2, 2, 6, 2, 2}),	// Si: [Ne] 3s2 3p2
    std::vector<int>({2, 2, 6, 2, 3}),	// P: [Ne] 3s2 3p3
    std::vector<int>({2, 2, 6, 2, 4}),	// S: [Ne] 3s2 3p4
    std::vector<int>({2, 2, 6, 2, 5}),	// Cl: [Ne] 3s2 3p5
    std::vector<int>({2, 2, 6, 2, 6}),	// Ar: [Ne] 3s2 3p6
    std::vector<int>({2, 2, 6, 2, 6, 1}),	// K: [Ar] 4s1
    std::vector<int>({2, 2, 6, 2, 6, 2}),	// Ca: [Ar] 4s2
    std::vector<int>({2, 2, 6, 2, 6, 2, 1}),	// Sc: [Ar] 4s2 3d1
    std::vector<int>({2, 2, 6, 2, 6, 2, 2}),	// Ti: [Ar] 4s2 3d2
    std::vector<int>({2, 2, 6, 2, 6, 2, 3}),	// V: [Ar] 4s2 3d3
    std::vector<int>({2, 2, 6, 2, 6, 1, 5}),	// Cr: [Ar] 4s1 3d5 *
    std::vector<int>({2, 2, 6, 2, 6, 2, 5}),	// Mn: [Ar] 4s2 3d5
    std::vector<int>({2, 2, 6, 2, 6, 2, 6}),	// Fe: [Ar] 4s2 3d6
    std::vector<int>({2, 2, 6, 2, 6, 2, 7}),	// Co: [Ar] 4s2 3d7
    std::vector<int>({2, 2, 6, 2, 6, 2, 8}),	// Ni: [Ar] 4s2 3d8
    std::vector<int>({2, 2, 6, 2, 6, 1, 10}),	// Cu: [Ar] 4s1 3d10 *
    std::vector<int>({2, 2, 6, 2, 6, 2, 10}),	// Zn: [Ar] 4s2 3d10
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 1}),	// Ga: [Ar] 4s2 3d10 4p1
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 2}),	// Ge: [Ar] 4s2 3d10 4p2
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 3}),	// As: [Ar] 4s2 3d10 4p3
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 4}),	// Se: [Ar] 4s2 3d10 4p4
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 5}),	// Br: [Ar] 4s2 3d10 4p5
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6}),	// Kr: [Ar] 4s2 3d10 4p6
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 1}),	// Rb: [Kr] 5s1
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2}),	// Sr: [Kr] 5s2
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 1}),	// Y: [Kr] 5s2 4d1
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 2}),	// Zr: [Kr] 5s2 4d2
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 3}),	// Nb: [Kr] 5s2 4d3
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 1, 5}),	// Mo: [Kr] 5s1 4d5 *
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 5}),	// Tc: [Kr] 5s2 4d5
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 6}),	// Ru: [Kr] 5s2 4d6
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 7}),	// Rh: [Kr] 5s2 4d7
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 0, 10}),	// Pd: [Kr] 4d10 *
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 1, 10}),	// Ag: [Kr] 5s1 4d10 *
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10}),	// Cd: [Kr] 5s2 4d10
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 1}),	// In: [Kr] 5s2 4d10 5p1
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 2}),	// Sn: [Kr] 5s2 4d10 5p2
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 3}),	// Sb: [Kr] 5s2 4d10 5p3
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 4}),	// Te: [Kr] 5s2 4d10 5p4
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 5}),	// I: [Kr] 5s2 4d10 5p5
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6}),	// Xe: [Kr] 5s2 4d10 5p6
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 1}),	// Cs: [Xe] 6s1
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2}),	// Ba: [Xe] 6s2
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 0, 1}),	// La: [Xe] 6s2 5d1 *
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 1, 1}),	// Ce: [Xe] 6s2 4f1 5d1 *
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 3}),	// Pr: [Xe] 6s2 4f3
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 4}),	// Nd: [Xe] 6s2 4f4
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 5}),	// Pm: [Xe] 6s2 4f5
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 6}),	// Sm: [Xe] 6s2 4f6
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 7}),	// Eu: [Xe] 6s2 4f7
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 7, 1}),	// Gd: [Xe] 6s2 4f7 5d1
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 9}),	// Tb: [Xe] 6s2 4f9
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 10}),	// Dy: [Xe] 6s2 4f10
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 11}),	// Ho: [Xe] 6s2 4f11
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 12}),	// Er: [Xe] 6s2 4f12
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 13}),	// Tm: [Xe] 6s2 4f13
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14}),	// Yb: [Xe] 6s2 4f14
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 1}),	// Lu: [Xe] 6s2 4f14 5d1
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 2}),	// Hf: [Xe] 6s2 4f14 5d2
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 3}),	// Ta: [Xe] 6s2 4f14 5d3
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 4}),	// W: [Xe] 6s2 4f14 5d4
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 5}),	// Re: [Xe] 6s2 4f14 5d5
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 6}),	// Os: [Xe] 6s2 4f14 5d6
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 7}),	// Ir: [Xe] 6s2 4f14 5d7
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 1, 14, 9}),	// Pt: [Xe] 6s1 4f14 5d9 *
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 1, 14, 10}),	// Au: [Xe] 6s1 4f14 5d10 *
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10}),	// Hg: [Xe] 6s2 4f14 5d10
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 1}),	// Tl: [Xe] 6s2 4f14 5d10 6p1
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 2}),	// Pb: [Xe] 6s2 4f14 5d10 6p2
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 3}),	// Bi: [Xe] 6s2 4f14 5d10 6p3
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 4}),	// Po: [Xe] 6s2 4f14 5d10 6p4
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 5}),	// At: [Xe] 6s2 4f14 5d10 6p5
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 6}),	// Rn: [Xe] 6s2 4f14 5d10 6p6
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 6, 1}),	// Fr: [Rn] 7s1
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 6, 2}),	// Ra: [Rn] 7s2
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 6, 2, 0, 1}),	// Ac: [Rn] 7s2 6d1 *
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 6, 2, 0, 2}),	// Th: [Rn] 7s2 6d2 *
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 6, 2, 2, 1}),	// Pa: [Rn] 7s2 5f2 6d1 *
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 6, 2, 3, 1}),	// U: [Rn] 7s2 5f3 6d1 *
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 6, 2, 4, 1}),	// Np: [Rn] 7s2 5f4 6d1 *
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 6, 2, 6}),	// Pu: [Rn] 7s2 5f6
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 6, 2, 7}),	// Am: [Rn] 7s2 5f7
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 6, 2, 7, 1}),	// Cm: [Rn] 7s2 5f7 6d1 *
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 6, 2, 9}),	// Bk: [Rn] 7s2 5f9
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 6, 2, 10}),	// Cf: [Rn] 7s2 5f10
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 6, 2, 11}),	// Es: [Rn] 7s2 5f11
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 6, 2, 12}),	// Fm: [Rn] 7s2 5f12
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 6, 2, 13}),	// Md: [Rn] 7s2 5f13
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 6, 2, 14}),	// No: [Rn] 7s2 5f14
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 6, 2, 14, 0, 1}),	// Lr: [Rn] 7s2 5f14 7p1 *
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 6, 2, 14, 2}),	// Rf: [Rn] 7s2 5f14 6d2
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 6, 2, 14, 3}),	// Db: [Rn] 7s2 5f14 6d3
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 6, 2, 14, 4}),	// Sg: [Rn] 7s2 5f14 6d4
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 6, 2, 14, 5}),	// Bh: [Rn] 7s2 5f14 6d5
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 6, 2, 14, 6}),	// Hs: [Rn] 7s2 5f14 6d6
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 6, 2, 14, 7}),	// Mt: [Rn] 7s2 5f14 6d7
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 6, 2, 14, 8}),	// Ds: [Rn] 7s2 5f14 6d8
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 6, 2, 14, 9}),	// Rg: [Rn] 7s2 5f14 6d9
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 6, 2, 14, 10}),	// Cn: [Rn] 7s2 5f14 6d10
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 6, 2, 14, 10, 1}),	// Nh: [Rn] 7s2 5f14 6d10 7p1
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 6, 2, 14, 10, 2}),	// Fl: [Rn] 7s2 5f14 6d10 7p2
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 6, 2, 14, 10, 3}),	// Mc: [Rn] 7s2 5f14 6d10 7p3
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 6, 2, 14, 10, 4}),	// Lv: [Rn] 7s2 5f14 6d10 7p4
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 6, 2, 14, 10, 5}),	// Tn: [Rn] 7s2 5f14 6d10 7p5
    std::vector<int>({2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 6, 2, 14, 10, 6})	// Og: [Rn] 7s2 5f14 6d10 7p6
  };

  if(Z < 1 || Z > 118) {
    throw new std::runtime_error("Only accepts atomic numbers between 1 and 118.");
  }
    
  GSConfig *out = new GSConfig(Z, gsconfigs[Z - 1]);
  return out;
}
