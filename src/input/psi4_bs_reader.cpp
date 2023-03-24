#include "psi4_bs_reader.hpp"
#include "../basis_sets/basis_set.hpp"
#include "../util/atom.hpp"
#include <stdexcept>
#include <cctype>
#include <cstring>

using namespace compchem;
using namespace std;

#define BUFF_SIZE 1024

std::vector<BasisOrbital *> *compchem::readPsi4file(std::FILE *fp, int Z, int charge) {

  std::rewind(fp); 
  std::vector<BasisOrbital *> *out = new std::vector<BasisOrbital *>();
  char buffer[BUFF_SIZE + 1];
  bool comment = false;

  int currz = 0, currl = 0, currcharge = 0, ngauss = 0, currgauss = 0;
  int end = 0;
  std::vector<double> alphas = std::vector<double>(),
    coefs = std::vector<double>(),
    coefs2 = std::vector<double>();
  double alpha, coef, coef2, norm;
  char currchar;
  int linenum = 1;

  enum {
    FIND_DATA,
    SYMBOL,
    CHARGE,
    ORBITAL_TYPE,
    GAUSSIANS,
    NORMALIZATION,
    PAIRS_FIRST,
    PAIRS_SECOND,
    TRIPLES_FIRST,
    TRIPLES_SECOND,
    TRIPLES_THIRD,
  } state;
  state = FIND_DATA;

  while(!feof(fp)) {
    
    currchar = fgetc(fp);

    // Check if the end of a comment.
    if(comment && currchar == '\n') {
      comment = false;
      end = 0;
      linenum++;
      continue;
    }

    // Check if in a comment but not at the end.
    if(comment) {
      end = 0;
      continue;
    }

    // End of a token.
    if(isspace(currchar) || (int) currchar == 0 || currchar == '!') {
      // Didn't get any input on the line, or the line is a comment.
      if(end == 0 && currchar == '!') {
	comment = true;
	continue;
      } else if(end == 0) {
	if(currchar == '\n') {
	  linenum++;
	}
	continue;
      }
      
      switch(state) {
      case FIND_DATA :
	if(buffer[0] == '*') {
	  state = SYMBOL;
	}
	break;
      case SYMBOL :
	currz = getZFromSymb(buffer);
	if(currz > 0) {
	  state = CHARGE;
	}
	break;
      case CHARGE :
	currcharge = stoi(buffer);
	state = ORBITAL_TYPE;
	break;
      case ORBITAL_TYPE :
	// Looking for an orbital, found next atom.
	if(buffer[0] == '*') {
	  if(currz == Z && currcharge == charge) {
	    return out;
	  } else {
	    // Look for the next symbol. Reset all parameters.
	    state = SYMBOL;
	    currz = 0;
	    currl = 0;
	    currcharge = 0;
	    ngauss = 0;
	    currgauss = 0;
	    break;
	  }
	}
	// Found an orbital.
	if(strncmp(buffer, "SP", 3) == 0) {
	  currl = -1;
	} else {
	  switch(buffer[0]) {
	  case 'S':
	    currl = 0;
	    break;
	  case 'P':
	    currl = 1;
	    break;
	  case 'D':
	    currl = 2;
	    break;
	  case 'F':
	    currl = 3;
	    break;
	  default :
	    if(isupper(buffer[0])) {
	      currl = 3 + buffer[0] - 'F';
	    } else {
	      for(auto b : *out) {
		delete b;
	      }
	      out->clear();
	      delete out;
	      throw new std::runtime_error("Could not recognize orbital symbol at line " + to_string(linenum) +  "\n");
	    }
	  }
	}

	coefs.clear();
	alphas.clear();
	coefs2.clear();
	currgauss = 0;
	
	state = GAUSSIANS;
	break;
      case GAUSSIANS :
	ngauss = stoi(buffer);
        state = NORMALIZATION;
	break;
      case NORMALIZATION :
	norm = stod(buffer);
	if(currl == -1) {
	  state = TRIPLES_FIRST;
	} else {
	  state = PAIRS_FIRST;
	}
	if(norm == 0) {
	  norm = 1;
	}
	currgauss = 0;
	break;
      case PAIRS_FIRST :
	alphas.push_back(stod(buffer));
	state = PAIRS_SECOND;
	break;
      case PAIRS_SECOND :
	coefs.push_back(stod(buffer) / norm);
	currgauss++;
	// Finished gathering orbital. Add it and all its subshells.
	if(currgauss >= ngauss) {
	  if(currz == Z && currcharge == charge) {
	    for(int ml = -currl; ml <= currl; ml++) {
	      out->push_back(new GaussianOrbital(currl, ml, coefs, alphas));
	    }
	  }
	  coefs.clear();
	  alphas.clear();
	  state = ORBITAL_TYPE;
	} else {
	  state = PAIRS_FIRST;
	}
	break;
      case TRIPLES_FIRST :
	alphas.push_back(stod(buffer));
	state = TRIPLES_SECOND;
	break;
      case TRIPLES_SECOND :
        coefs.push_back(stod(buffer) / norm);
	state = TRIPLES_THIRD;
	break;
      case TRIPLES_THIRD :
	coefs2.push_back(stod(buffer) / norm);
	currgauss++;
	// Finish gathering orbitals. Add the s and p orbitals.
	if(currgauss >= ngauss) {
	  if(currz == Z && currcharge == charge) {
	    out->push_back(new GaussianOrbital(0, 0, coefs, alphas));
	    for(int ml = -1; ml <= 1; ml++) {
	      out->push_back(new GaussianOrbital(1, ml, coefs2, alphas));
	    }
	  }
	  alphas.clear();
	  coefs.clear();
	  coefs2.clear();
	  state = ORBITAL_TYPE;
	} else {
	  state = TRIPLES_FIRST;
	}
	break;
      }
      if(currchar == '\n') {
	linenum++;
      }
      end = 0;
    } else {
      // Not end of a token. Add character.
      if(end == BUFF_SIZE) {
	for(auto b : *out) {
	  delete b;
	}
	out->clear();
	delete out;
	throw new std::runtime_error("Malformed line at line " + to_string(linenum) + "\n");
      }
      buffer[end] = currchar;
      buffer[end + 1] = 0;
      end++;
    }
  }

  for(auto b : *out) {
    delete b;
  }
  out->clear();
  delete out;
  throw new std::runtime_error("Basis set not found!\n");
}
