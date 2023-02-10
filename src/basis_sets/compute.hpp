/*
 * compute.hpp
 *
 *  Created on: Feb 7, 2023
 *      Author: connor
 */

#ifndef __COMPUTE_HPP__
#define __COMPUTE_HPP__

#include "../ints/gaussian.hpp"
#include <stdlib.h>
#include <stdio.h>

namespace compchem {

class ShellFactory {
private :
	static ShellFactory *singleton = nullptr;

	ShellFactory();
public :
	~ShellFactory() = default;

	static ShellFactory &getSingleton();

	GaussianSum &makeShell(int l, int m, double alpha);
};

class OrbitalComputer {
public :
	virtual ~OrbitalComputer() = 0;

	GaussianSum *compute(int Z, int *norbs);
};

class OrbitalPrecomputer {

public:
	virtual ~OrbitalPrecomputer() = 0;

	void compute_orbital(FILE *out, int *Zs, int nzs, OrbitalComputer &computer);
};

class OrbitalFactory {
private :
	static OrbitalFactory *singleton;

	OrbitalFactory();
public:
	static OrbitalFactory &getsingleton();

	// Read from files.
	void makePsi4(int Z, GaussianSum *orbitals, int *norbs, const std::array<double, 3> &center, FILE *fp);
	void makePrecomp(int Z, GaussianSum *orbitals, int *norbs, const std::array<double, 3> &center, FILE *fp);
	void makeDirect(int Z, GaussianSum *orbitals, int *norbs, const std::array<double, 3> &center, OrbitalComputer &computer);
};

int shells(int Z);

}

#endif


