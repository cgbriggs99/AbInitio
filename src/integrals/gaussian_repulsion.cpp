 
#include <extramath.h>
#include "../basis_sets/basis_set.hpp"
#include "../util/atom.hpp"
#include <array>
#include <cmath>
#include <vector>
#include <deque>
#include <map>
#include "../util/polynomial.hpp"
#include "integrals.hpp"

using namespace compchem;    

template<typename T>
static bool contains(const T &item, const std::deque<T> &col) {
  for(int i = 0; i < col.size(); i++) {
    if(col.at(i) == item) {
      return true;
    }
  }
  return false;
}

// Figure out which e0f0 entries need to be found.
static void sieve_abcd(const std::array<int, 12> &index,
		       std::deque<std::array<int, 6> > &e0f0) {
  std::deque<std::array<int, 12> > open, closed;
  open.push_back(index);
  
  while(open.size() != 0) {
    const std::array<int, 12> &curr_ind = open.front();
    // Check that the index is valid.
    for(int i = 0; i < 12; i++) {
      if(curr_ind[i] < 0) {
	// Finish out the loop.
	goto LOOP_1;
      }
    }

    // Don't redo any already done.
    if(contains(curr_ind, closed)) {
      open.pop_front();
      continue;
    }
    
    if(curr_ind[3] == 0) {
      if(curr_ind[4] == 0) {
	if(curr_ind[5] == 0) {
	  if(curr_ind[9] == 0) {
	    if(curr_ind[10] == 0) {
	      if(curr_ind[11] == 0) {
		// Index is an (e0|f0) index.
		std::array<int, 6> ind1 = {
		  curr_ind[0], curr_ind[1], curr_ind[2],
		  curr_ind[6], curr_ind[7], curr_ind[8]};
		if(!contains(ind1, e0f0)) {
		  e0f0.push_back(ind1);
		}
	      } else {
		std::array<int, 12> ind1 = {
		  curr_ind[0], curr_ind[1], curr_ind[2],
		  curr_ind[3], curr_ind[4], curr_ind[5],
		  curr_ind[6], curr_ind[7], curr_ind[8] + 1,
		  curr_ind[9], curr_ind[10], curr_ind[11] - 1},
		  ind2 = {
		    curr_ind[0], curr_ind[1], curr_ind[2],
		    curr_ind[3], curr_ind[4], curr_ind[5],
		    curr_ind[6], curr_ind[7], curr_ind[8],
		    curr_ind[9], curr_ind[10], curr_ind[11] - 1};
		open.push_back(ind1);
		open.push_back(ind2);
	      }
	    } else {
	      std::array<int, 12> ind1 = {
		curr_ind[0], curr_ind[1], curr_ind[2],
		curr_ind[3], curr_ind[4], curr_ind[5],
		curr_ind[6], curr_ind[7] + 1, curr_ind[8],
		curr_ind[9], curr_ind[10] - 1, curr_ind[11]},
		ind2 = {
		  curr_ind[0], curr_ind[1], curr_ind[2],
		  curr_ind[3], curr_ind[4], curr_ind[5],
		  curr_ind[6], curr_ind[7], curr_ind[8],
		  curr_ind[9], curr_ind[10] - 1, curr_ind[11]};
	      open.push_back(ind1);
	      open.push_back(ind2);
	    }
	  } else {
	    std::array<int, 12> ind1 = {
	      curr_ind[0], curr_ind[1], curr_ind[2],
	      curr_ind[3], curr_ind[4], curr_ind[5],
	      curr_ind[6] + 1, curr_ind[7], curr_ind[8],
	      curr_ind[9] - 1, curr_ind[10], curr_ind[11]},
	      ind2 = {
		curr_ind[0], curr_ind[1], curr_ind[2],
		curr_ind[3], curr_ind[4], curr_ind[5],
		curr_ind[6], curr_ind[7], curr_ind[8],
		curr_ind[9] - 1, curr_ind[10], curr_ind[11]};
	      open.push_back(ind1);
	      open.push_back(ind2);
	  }
	} else {
	  std::array<int, 12> ind1 = {
	    curr_ind[0], curr_ind[1], curr_ind[2] + 1,
	    curr_ind[3], curr_ind[4], curr_ind[5] - 1,
	    curr_ind[6], curr_ind[7], curr_ind[8],
	    curr_ind[9], curr_ind[10], curr_ind[11]},
	    ind2 = {
	      curr_ind[0], curr_ind[1], curr_ind[2],
	      curr_ind[3], curr_ind[4], curr_ind[5] - 1,
	      curr_ind[6], curr_ind[7], curr_ind[8],
	      curr_ind[9], curr_ind[10], curr_ind[11]};
	  open.push_back(ind1);
	  open.push_back(ind2);
	}
      } else {
	std::array<int, 12> ind1 = {
	  curr_ind[0], curr_ind[1] + 1, curr_ind[2],
	  curr_ind[3], curr_ind[4] - 1, curr_ind[5],
	  curr_ind[6], curr_ind[7], curr_ind[8],
	  curr_ind[9], curr_ind[10], curr_ind[11]},
	  ind2 = {
	    curr_ind[0], curr_ind[1], curr_ind[2],
	    curr_ind[3], curr_ind[4] - 1, curr_ind[5],
	    curr_ind[6], curr_ind[7], curr_ind[8],
	    curr_ind[9], curr_ind[10], curr_ind[11]};
	open.push_back(ind1);
	open.push_back(ind2);
      }
    } else {
      std::array<int, 12> ind1 = {
	curr_ind[0] + 1, curr_ind[1], curr_ind[2],
	curr_ind[3] - 1, curr_ind[4], curr_ind[5],
	curr_ind[6], curr_ind[7], curr_ind[8],
	curr_ind[9], curr_ind[10], curr_ind[11]},
	ind2 = {
	  curr_ind[0], curr_ind[1], curr_ind[2],
	  curr_ind[3] - 1, curr_ind[4], curr_ind[5],
	  curr_ind[6], curr_ind[7], curr_ind[8],
	  curr_ind[9], curr_ind[10], curr_ind[11]};
      open.push_back(ind1);
      open.push_back(ind2);
    }

  LOOP_1:;
    closed.push_back(curr_ind);
    open.pop_front();
  }
}

// Find the e0quv indices.
static void sieve_e0f0(const std::deque<std::array<int, 6> > &in,
		       std::deque<std::array<int, 8> > &e0quv,
		       std::deque<std::array<int, 6> > &e0q) {
  std::deque<std::array<int, 11> > open, closed;
  for(int i = 0; i < in.size(); i++) {
    const std::array<int, 6> &curr_ind = in.at(i);
    std::array<int, 11> ind = {
      curr_ind[0], curr_ind[1], curr_ind[2],
      curr_ind[3], curr_ind[4], curr_ind[5],
      0, 0, 0, 0, 0};
    open.push_back(ind);
  }

  while(open.size() != 0) {
    const std::array<int, 11> &curr_ind = open.front();

    // Check that the index is valid.
    for(int i = 0; i < 11; i++) {
      if(curr_ind[i] < 0) {
	// If index is not valid, finish out the loop.
	goto LOOP_2;
      }
    }

    // Don't redo any that have already been done.
    if(contains(curr_ind, closed)) {
      open.pop_front();
      continue;
    }

    if(curr_ind[3] == 0) {
      if(curr_ind[4] == 0) {
	if(curr_ind[5] == 0) {
	  std::array<int, 8> ind1 = {
	    curr_ind[0], curr_ind[1], curr_ind[2],
	    curr_ind[6], curr_ind[7], curr_ind[8],
	    curr_ind[9], curr_ind[10]};
	  if(!contains(ind1, e0quv)) {
	    e0quv.push_back(ind1);
	  }
	  std::array<int, 6> ind2 = {
	    curr_ind[0], curr_ind[1], curr_ind[2],
	    curr_ind[6], curr_ind[7], curr_ind[8]};
	  if(!contains(ind2, e0q)) {
	    e0q.push_back(ind2);
	  }
	} else {
	  std::array<int, 11> ind1 = {
	    curr_ind[0], curr_ind[1], curr_ind[2],
	    curr_ind[3], curr_ind[4], curr_ind[5] - 1,
	    curr_ind[6], curr_ind[7], curr_ind[8] - 1,
	    curr_ind[9], curr_ind[10]},
	    ind2 = {
	    curr_ind[0], curr_ind[1], curr_ind[2],
	    curr_ind[3], curr_ind[4], curr_ind[5] - 1,
	    curr_ind[6], curr_ind[7], curr_ind[8],
	    curr_ind[9] + 1, curr_ind[10] + 1},
	    ind3 = {
	    curr_ind[0], curr_ind[1], curr_ind[2],
	    curr_ind[3], curr_ind[4], curr_ind[5] - 1,
	    curr_ind[6], curr_ind[7], curr_ind[8] + 1,
	    curr_ind[9], curr_ind[10] + 1};
	  open.push_back(ind1);
	  open.push_back(ind2);
	  open.push_back(ind3);
	}
      } else {
	std::array<int, 11> ind1 = {
	  curr_ind[0], curr_ind[1], curr_ind[2],
	  curr_ind[3], curr_ind[4] - 1, curr_ind[5],
	  curr_ind[6], curr_ind[7] - 1, curr_ind[8],
	  curr_ind[9], curr_ind[10]},
	  ind2 = {
	    curr_ind[0], curr_ind[1], curr_ind[2],
	    curr_ind[3], curr_ind[4] - 1, curr_ind[5],
	    curr_ind[6], curr_ind[7], curr_ind[8],
	    curr_ind[9] + 1, curr_ind[10] + 1},
	  ind3 = {
	    curr_ind[0], curr_ind[1], curr_ind[2],
	    curr_ind[3], curr_ind[4] - 1, curr_ind[5],
	    curr_ind[6], curr_ind[7] + 1, curr_ind[8],
	    curr_ind[9], curr_ind[10] + 1};
	  open.push_back(ind1);
	  open.push_back(ind2);
	  open.push_back(ind3);
      }
    } else {
      std::array<int, 11> ind1 = {
	curr_ind[0], curr_ind[1], curr_ind[2],
	curr_ind[3] - 1, curr_ind[4], curr_ind[5],
	curr_ind[6] - 1, curr_ind[7], curr_ind[8],
	curr_ind[9], curr_ind[10]},
	ind2 = {
	  curr_ind[0], curr_ind[1], curr_ind[2],
	  curr_ind[3] - 1, curr_ind[4], curr_ind[5],
	  curr_ind[6], curr_ind[7], curr_ind[8],
	  curr_ind[9] + 1, curr_ind[10] + 1},
	ind3 = {
	  curr_ind[0], curr_ind[1], curr_ind[2],
	  curr_ind[3] - 1, curr_ind[4], curr_ind[5],
	  curr_ind[6] + 1, curr_ind[7], curr_ind[8],
	  curr_ind[9], curr_ind[10] + 1};
      open.push_back(ind1);
      open.push_back(ind2);
      open.push_back(ind3);
    }
  LOOP_2:;
    closed.push_back(curr_ind);
    open.pop_front();
  }
}

// Compute [r]_uv indices.
static void sieve_pquv(const std::deque<std::array<int, 6> > &in,
		       std::deque<std::array<int, 5> > &ruv,
		       std::deque<std::array<int, 3> > &r) {
  std::deque<std::array<int, 11> > open, closed;

  for(int i = 0; i < in.size(); i++) {
    const std::array<int, 6> &curr_ind = in.at(i);
    std::array<int, 11> ind = {
      curr_ind[0], curr_ind[1], curr_ind[2],
      0, 0, 0,
      curr_ind[3], curr_ind[4], curr_ind[5],
      0, 0};
    open.push_back(ind);
  }

  while(open.size() != 0) {
    const std::array<int, 11> &curr_ind = open.front();

    // Check that the index is valid.
    for(int i = 0; i < 11; i++) {
      if(curr_ind[i] < 0) {
	goto LOOP_3;
      }
    }

    if(contains(curr_ind, closed)) {
      open.pop_front();
      continue;
    }

    if(curr_ind[0] == 0) {
      if(curr_ind[1] == 0) {
	if(curr_ind[2] == 0) {
	  std::array<int, 5> ind1 = {
	    curr_ind[3] + curr_ind[6],
	    curr_ind[4] + curr_ind[7],
	    curr_ind[5] + curr_ind[8],
	    curr_ind[9], curr_ind[10]};
	  if(!contains(ind1, ruv)) {
	    ruv.push_back(ind1);
	  }
	  std::array<int, 3> ind2 = {
	    curr_ind[3] + curr_ind[6],
	    curr_ind[4] + curr_ind[7],
	    curr_ind[5] + curr_ind[8]};
	  if(!contains(ind2, r)) {
	    r.push_back(ind2);
	  }
	} else {
	  std::array<int, 11> ind1 = {
	    curr_ind[0], curr_ind[1], curr_ind[2] - 1,
	    curr_ind[3], curr_ind[4], curr_ind[5] - 1,
	    curr_ind[6], curr_ind[7], curr_ind[8],
	    curr_ind[9], curr_ind[10]},
	    ind2 = {
	    curr_ind[0], curr_ind[1], curr_ind[2] - 1,
	    curr_ind[3], curr_ind[4], curr_ind[5],
	    curr_ind[6], curr_ind[7], curr_ind[8],
	    curr_ind[9] + 1, curr_ind[10] + 1},
	    ind3 = {
	    curr_ind[0], curr_ind[1], curr_ind[2] - 1,
	    curr_ind[3], curr_ind[4], curr_ind[5] + 1,
	    curr_ind[6], curr_ind[7], curr_ind[8],
	    curr_ind[9], curr_ind[10] + 1};
	  open.push_back(ind1);
	  open.push_back(ind2);
	  open.push_back(ind3);
	}
      } else {
	std::array<int, 11> ind1 = {
	  curr_ind[0], curr_ind[1] - 1, curr_ind[2],
	  curr_ind[3], curr_ind[4] - 1, curr_ind[5],
	  curr_ind[6], curr_ind[7], curr_ind[8],
	  curr_ind[9], curr_ind[10]},
	  ind2 = {
	    curr_ind[0], curr_ind[1] - 1, curr_ind[2],
	    curr_ind[3], curr_ind[4], curr_ind[5],
	    curr_ind[6], curr_ind[7], curr_ind[8],
	    curr_ind[9] + 1, curr_ind[10] + 1},
	  ind3 = {
	    curr_ind[0], curr_ind[1] - 1, curr_ind[2],
	    curr_ind[3], curr_ind[4] + 1, curr_ind[5],
	    curr_ind[6], curr_ind[7], curr_ind[8],
	    curr_ind[9], curr_ind[10] + 1};
	open.push_back(ind1);
	open.push_back(ind2);
	open.push_back(ind3);
      }
    } else {
      std::array<int, 11> ind1 = {
	curr_ind[0] - 1, curr_ind[1], curr_ind[2],
	curr_ind[3] - 1, curr_ind[4], curr_ind[5],
	curr_ind[6], curr_ind[7], curr_ind[8],
	curr_ind[9], curr_ind[10]},
	ind2 = {
	  curr_ind[0] - 1, curr_ind[1], curr_ind[2],
	  curr_ind[3], curr_ind[4], curr_ind[5],
	  curr_ind[6], curr_ind[7], curr_ind[8],
	  curr_ind[9] + 1, curr_ind[10] + 1},
	ind3 = {
	  curr_ind[0] - 1, curr_ind[1], curr_ind[2],
	  curr_ind[3] + 1, curr_ind[4], curr_ind[5],
	  curr_ind[6], curr_ind[7], curr_ind[8],
	  curr_ind[9], curr_ind[10] + 1};
      open.push_back(ind1);
      open.push_back(ind2);
      open.push_back(ind3);
    }
  LOOP_3:;
    closed.push_back(curr_ind);
    open.pop_front();
  }
}

static void sieve_0m(const std::deque<std::array<int, 3> > &in,
		     std::deque<int> &out) {
  std::deque<std::array<int, 4> > open, closed;
  for(int i = 0; i < in.size(); i++) {
    const std::array<int, 3> &curr_ind = in.at(i);

    std::array<int, 4> ind = {
      curr_ind[0], curr_ind[1], curr_ind[2], 0};
    open.push_back(ind);
  }

  while(open.size() != 0) {
    const std::array<int, 4> &curr_ind = open.front();

    for(int i = 0; i < 4; i++) {
      if(curr_ind[i] < 0) {
	goto LOOP_4;
      }
    }

    if(contains(curr_ind, closed)) {
      open.pop_front();
      continue;
    }

    if(curr_ind[0] == 0) {
      if(curr_ind[1] == 0) {
	if(curr_ind[2] == 0) {
	  if(!contains(curr_ind[3], out)) {
	    out.push_back(curr_ind[3]);
	  }
	} else {
	  std::array<int, 4> ind1 = {
	    curr_ind[0], curr_ind[1], curr_ind[2] - 1, curr_ind[3] + 1},
	    ind2 = {
	      curr_ind[0], curr_ind[1], curr_ind[2] - 2, curr_ind[3] + 1};
	  open.push_back(ind1);
	  open.push_back(ind2);
	}
      } else {
	std::array<int, 4> ind1 = {
	  curr_ind[0], curr_ind[1] - 1, curr_ind[2], curr_ind[3] + 1},
	  ind2 = {
	    curr_ind[0], curr_ind[1] - 2, curr_ind[2], curr_ind[3] + 1};
	open.push_back(ind1);
	open.push_back(ind2);
      }
    } else {
      std::array<int, 4> ind1 = {
	curr_ind[0] - 1, curr_ind[1], curr_ind[2], curr_ind[3] + 1},
	ind2 = {
	  curr_ind[0] - 2, curr_ind[1], curr_ind[2], curr_ind[3] + 1};
      open.push_back(ind1);
      open.push_back(ind2);
    }
  LOOP_4:;
    closed.push_back(curr_ind);
    open.pop_front();
  }
}

static double compute_rm(const std::array<int, 4> &index,
			 std::map<std::array<int, 4>, double> &ints,
			 const std::map<int, double> &m_ints,
			 double Rx, double Ry, double Rz) {
  if(ints.count(index) != 0) {
    return ints.at(index);
  }

  for(int i = 0; i < 4; i++) {
    if(index[i] < 0) {
      return 0;
    }
  }

  if(index[0] == 0) {
    if(index[1] == 0) {
      if(index[2] == 0) {
	return m_ints.at(index[3]);
      } else {
	std::array<int, 4> ind1 = {
	  index[0], index[1], index[2] - 1, index[3] + 1},
	  ind2 = {
	    index[0], index[1], index[2] - 2, index[3] + 1};
	double res = Rz * compute_rm(ind1, ints, m_ints, Rx, Ry, Rz) -
	  (index[2] - 1) * compute_rm(ind2, ints, m_ints, Rx, Ry, Rz);
	ints[index] = res;
	return res;
      }
    } else {
      std::array<int, 4> ind1 = {
	index[0], index[1] - 1, index[2], index[3] + 1},
	ind2 = {
	  index[0], index[1] - 2, index[2], index[3] + 1};
      double res = Ry * compute_rm(ind1, ints, m_ints, Rx, Ry, Rz) -
	(index[1] - 1) * compute_rm(ind2, ints, m_ints, Rx, Ry, Rz);
      ints[index] = res;
      return res;
    }
  } else {
    std::array<int, 4> ind1 = {
      index[0] - 1, index[1], index[2], index[3] + 1},
      ind2 = {
	index[0] - 2, index[1], index[2], index[3] + 1};
    double res = Rx * compute_rm(ind1, ints, m_ints, Rx, Ry, Rz) -
      (index[0] - 1) * compute_rm(ind2, ints, m_ints, Rx, Ry, Rz);
    ints[index] = res;
    return res;
  }
}

static double compute_apquv(const std::array<int, 11> &index,
			    std::map<std::array<int, 11>, double> &ints,
			    const std::map<std::array<int, 5>, double> &pquv,
			    const std::array<double, 3> &c1,
			    const std::array<double, 3> &c2) {
  if(ints.count(index) != 0) {
    return ints.at(index);
  }

  for(int i = 0; i < 11; i++) {
    if(index[i] < 0) {
      return 0;
    }
  }

  if(index[0] == 0) {
    if(index[1] == 0) {
      if(index[2] == 0) {
	int sign = ((index[6] + index[7] + index[8]) % 2) ? -1 : 1;
	std::array<int, 5> ind1 = {
	  index[3] + index[6],
	  index[4] + index[7],
	  index[5] + index[8],
	  index[9], index[10]};
        return sign * pquv.at(ind1);
      } else {
	std::array<int, 11> ind1 = {
	  index[0], index[1], index[2] - 1,
	  index[3], index[4], index[5] - 1,
	  index[6], index[7], index[8],
	  index[9], index[10]},
	  ind2 = {
	  index[0], index[1], index[2] - 1,
	  index[3], index[4], index[5],
	  index[6], index[7], index[8],
	  index[9] + 1, index[10] + 1},
	  ind3 = {
	  index[0], index[1], index[2] - 1,
	  index[3], index[4], index[5] + 1,
	  index[6], index[7], index[8],
	  index[9], index[10] + 1};
	double ret = index[5] * compute_apquv(ind1, ints, pquv, c1, c2) -
	  (c1[2] - c2[2]) * compute_apquv(ind2, ints, pquv, c1, c2) +
	  compute_apquv(ind3, ints, pquv, c1, c2);
	ints[index] = ret;
	return ret;
      }
    } else {
      std::array<int, 11> ind1 = {
	index[0], index[1] - 1, index[2],
	index[3], index[4] - 1, index[5],
	index[6], index[7], index[8],
	index[9], index[10]},
	ind2 = {
	  index[0], index[1] - 1, index[2],
	  index[3], index[4], index[5],
	  index[6], index[7], index[8],
	  index[9] + 1, index[10] + 1},
	ind3 = {
	  index[0], index[1] - 1, index[2],
	  index[3], index[4] + 1, index[5],
	  index[6], index[7], index[8],
	  index[9], index[10] + 1};
      double ret = index[4] * compute_apquv(ind1, ints, pquv, c1, c2) -
	(c1[1] - c2[1]) * compute_apquv(ind2, ints, pquv, c1, c2) +
	compute_apquv(ind3, ints, pquv, c1, c2);
      ints[index] = ret;
      return ret;
    }
  } else {
    std::array<int, 11> ind1 = {
      index[0] - 1, index[1], index[2],
      index[3] - 1, index[4], index[5],
      index[6], index[7], index[8],
      index[9], index[10]},
      ind2 = {
	index[0] - 1, index[1], index[2],
	index[3], index[4], index[5],
	index[6], index[7], index[8],
	index[9] + 1, index[10] + 1},
      ind3 = {
	index[0] - 1, index[1], index[2],
	index[3] + 1, index[4], index[5],
	index[6], index[7], index[8],
	index[9], index[10] + 1};
    double ret = index[3] * compute_apquv(ind1, ints, pquv, c1, c2) -
      (c1[0] - c2[0]) * compute_apquv(ind2, ints, pquv, c1, c2) +
      compute_apquv(ind3, ints, pquv, c1, c2);
    ints[index] = ret;
    return ret;
  }
}

static double compute_e0cquv(const std::array<int, 11> &index,
			     std::map<std::array<int, 11>, double> &ints,
			     const std::map<std::array<int, 8>, double> &e0quv,
			     const std::array<double, 3> &c3,
			     const std::array<double, 3> &c4) {
  if(ints.count(index) != 0) {
    return ints.at(index);
  }

  for(int i = 0; i < 11; i++) {
    if(index[i] < 0) {
      return 0;
    }
  }

  if(index[3] == 0) {
    if(index[4] == 0) {
      if(index[5] == 0) {
	std::array<int, 8> ind1 = {
	  index[0], index[1], index[2],
	  index[6], index[7], index[8],
	  index[9], index[10]};
	return e0quv.at(ind1);
      } else {
	std::array<int, 11> ind1 = {
	  index[0], index[1], index[2],
	  index[3], index[4], index[5] - 1,
	  index[6], index[7], index[8] - 1,
	  index[9], index[10]},
	  ind2 = {
	    index[0], index[1], index[2],
	    index[3], index[4], index[5] - 1,
	    index[6], index[7], index[8],
	    index[9] + 1, index[10] + 1},
	  ind3 = {
	    index[0], index[1], index[2],
	    index[3], index[4], index[5] - 1,
	    index[6], index[7], index[8] + 1,
	    index[9], index[10] + 1};
	double ret = index[8] * compute_e0cquv(ind1, ints, e0quv, c3, c4) -
	  (c3[2] - c4[2]) * compute_e0cquv(ind2, ints, e0quv, c3, c4) +
	  compute_e0cquv(ind3, ints, e0quv, c3, c4);
	ints[index] = ret;
	return ret;
      }
    } else {
      std::array<int, 11> ind1 = {
	index[0], index[1], index[2],
	index[3], index[4] - 1, index[5],
	index[6], index[7] - 1, index[8],
	index[9], index[10]},
	ind2 = {
	  index[0], index[1], index[2],
	  index[3], index[4] - 1, index[5],
	  index[6], index[7], index[8],
	  index[9] + 1, index[10] + 1},
	ind3 = {
	  index[0], index[1], index[2],
	  index[3], index[4] - 1, index[5],
	  index[6], index[7] + 1, index[8],
	  index[9], index[10] + 1};
      double ret = index[7] * compute_e0cquv(ind1, ints, e0quv, c3, c4) -
	(c3[1] - c4[1]) * compute_e0cquv(ind2, ints, e0quv, c3, c4) +
	compute_e0cquv(ind3, ints, e0quv, c3, c4);
      ints[index] = ret;
      return ret;
    }
  } else {
    std::array<int, 11> ind1 = {
      index[0], index[1], index[2],
      index[3] - 1, index[4], index[5],
      index[6] - 1, index[7], index[8],
      index[9], index[10]},
      ind2 = {
	index[0], index[1], index[2],
	index[3] - 1, index[4], index[5],
	index[6], index[7], index[8],
	index[9] + 1, index[10] + 1},
      ind3 = {
	index[0], index[1], index[2],
	index[3] - 1, index[4], index[5],
	index[6] + 1, index[7], index[8],
	index[9], index[10] + 1};
    double ret = index[6] * compute_e0cquv(ind1, ints, e0quv, c3, c4) -
      (c3[0] - c4[0]) * compute_e0cquv(ind2, ints, e0quv, c3, c4) +
      compute_e0cquv(ind3, ints, e0quv, c3, c4);
    ints[index] = ret;
    return ret;
  }
}

static double compute_abcd(const std::array<int, 12> &index,
			   std::map<std::array<int, 12>, double> &ints,
			   const std::map<std::array<int, 6>, double> &e0f0,
			   const std::array<double, 3> &c1,
			   const std::array<double, 3> &c2,
			   const std::array<double, 3> &c3,
			   const std::array<double, 3> &c4) {
  if(ints.count(index) != 0) {
    return ints.at(index);
  }

  for(int i = 0; i < 12; i++) {
    if(index[i] < 0) {
      return 0;
    }
  }

  if(index[3] == 0) {
    if(index[4] == 0) {
      if(index[5] == 0) {
	if(index[9] == 0) {
	  if(index[10] == 0) {
	    if(index[11] == 0) {
	      std::array<int, 6> ind1 = {
		index[0], index[1], index[2],
		index[6], index[7], index[8]};
	      return e0f0.at(ind1);
	    } else {
	      std::array<int, 12> ind1 = {
		index[0], index[1], index[2],
		index[3], index[4], index[5],
		index[6], index[7], index[8] + 1,
		index[9], index[10], index[11] - 1},
		ind2 = {
		  index[0], index[1], index[2],
		  index[3], index[4], index[5],
		  index[6], index[7], index[8],
		  index[9], index[10], index[11] - 1};
	      double ret = compute_abcd(ind1, ints, e0f0, c1, c2, c3, c4) +
		(c3[2] - c4[2]) * compute_abcd(ind2, ints, e0f0, c1, c2, c3, c4);
	      ints[index] = ret;
	      return ret;
	    }
	  } else {
	    std::array<int, 12> ind1 = {
	      index[0], index[1], index[2],
	      index[3], index[4], index[5],
	      index[6], index[7] + 1, index[8],
	      index[9], index[10] - 1, index[11]},
	      ind2 = {
		index[0], index[1], index[2],
		index[3], index[4], index[5],
		index[6], index[7], index[8],
		index[9], index[10] - 1, index[11]};
	    double ret = compute_abcd(ind1, ints, e0f0, c1, c2, c3, c4) +
	      (c3[1] - c4[1]) * compute_abcd(ind2, ints, e0f0, c1, c2, c3, c4);
	    ints[index] = ret;
	    return ret;
	  }
	} else {
	  std::array<int, 12> ind1 = {
	    index[0], index[1], index[2],
	    index[3], index[4], index[5],
	    index[6] + 1, index[7], index[8],
	    index[9] - 1, index[10], index[11]},
	    ind2 = {
	      index[0], index[1], index[2],
	      index[3], index[4], index[5],
	      index[6], index[7], index[8],
	      index[9] - 1, index[10], index[11]};
	  double ret = compute_abcd(ind1, ints, e0f0, c1, c2, c3, c4) +
	    (c3[0] - c4[0]) * compute_abcd(ind2, ints, e0f0, c1, c2, c3, c4);
	  ints[index] = ret;
	  return ret;
	}
      } else {
	std::array<int, 12> ind1 = {
	  index[0], index[1], index[2] + 1,
	  index[3], index[4], index[5] - 1,
	  index[6], index[7], index[8],
	  index[9], index[10], index[11]},
	  ind2 = {
	    index[0], index[1], index[2],
	    index[3], index[4], index[5] - 1,
	    index[6], index[7], index[8],
	    index[9], index[10], index[11]};
	double ret = compute_abcd(ind1, ints, e0f0, c1, c2, c3, c4) +
	  (c1[2] - c2[2]) * compute_abcd(ind2, ints, e0f0, c1, c2, c3, c4);
	ints[index] = ret;
	return ret;
      }
    } else {
      std::array<int, 12> ind1 = {
	index[0], index[1] + 1, index[2],
	index[3], index[4] - 1, index[5],
	index[6], index[7], index[8],
	index[9], index[10], index[11]},
	ind2 = {
	  index[0], index[1], index[2],
	  index[3], index[4] - 1, index[5],
	  index[6], index[7], index[8],
	  index[9], index[10], index[11]};
      double ret = compute_abcd(ind1, ints, e0f0, c1, c2, c3, c4) +
	(c1[1] - c2[1]) * compute_abcd(ind2, ints, e0f0, c1, c2, c3, c4);
      ints[index] = ret;
      return ret;
    }
  } else {
    std::array<int, 12> ind1 = {
      index[0] + 1, index[1], index[2],
      index[3] - 1, index[4], index[5],
      index[6], index[7], index[8],
      index[9], index[10], index[11]},
      ind2 = {
	index[0], index[1], index[2],
	index[3] - 1, index[4], index[5],
	index[6], index[7], index[8],
	index[9], index[10], index[11]};
    double ret = compute_abcd(ind1, ints, e0f0, c1, c2, c3, c4) +
      (c1[0] - c2[0]) * compute_abcd(ind2, ints, e0f0, c1, c2, c3, c4);
    ints[index] = ret;
    return ret;
  }
}
  
double AnalyticIntegral::rep_integral(const int *pows1, const int *pows2, const int *pows3,
				      const int *pows4,
				      const std::array<double, 3> &c1,
				      const std::array<double, 3> &c2,
				      const std::array<double, 3> &c3,
				      const std::array<double, 3> &c4,
				      const GaussianOrbital *o1,
				      const GaussianOrbital *o2,
				      const GaussianOrbital *o3,
				      const GaussianOrbital *o4) const {

  // Find which integrals need to be found.
  std::array<int, 12> index = {
    pows1[0], pows1[1], pows1[2],
    pows2[0], pows2[1], pows2[2],
    pows3[0], pows3[1], pows3[2],
    pows4[0], pows4[1], pows4[2]};

  std::deque<std::array<int, 6> > e0f0_inds, e0q_inds;
  std::deque<std::array<int, 8> > e0quv_inds;
  std::deque<std::array<int, 5> > pquv_inds;
  std::deque<std::array<int, 3> > r_inds;
  std::deque<int> m_inds;

  sieve_abcd(index, e0f0_inds);
  sieve_e0f0(e0f0_inds, e0quv_inds, e0q_inds);
  sieve_pquv(e0q_inds, pquv_inds, r_inds);
  sieve_0m(r_inds, m_inds);

  std::map<std::array<int, 8>, double> e0quv;
  
  for(int i = 0; i < o3->getnterms(); i++) {
    for(int j = 0; j < o4->getnterms(); j++) {
      double gamma = o3->getalpha(i), delta = o4->getalpha(j);
      double eta = gamma + delta;
      double Qx = (gamma * c3[0] + delta * c4[0]) / eta,
	Qy = (gamma * c3[1] + delta * c4[1]) / eta,
	Qz = (gamma * c3[2] + delta * c4[2]) / eta;
      double Kq = M_SQRT2 * std::pow(M_PI, 1.25) / eta *
	std::exp(-gamma * delta / eta *
		 ((c3[0] - c4[0]) * (c3[0] - c4[0]) +
		  (c3[1] - c4[1]) * (c3[1] - c4[1]) +
		  (c3[2] - c4[2]) * (c3[2] - c4[2])));
      double Dq = o3->getcoef(i) * o3->getnorm(i) *
	o4->getcoef(j) * o4->getnorm(j);

      std::map<std::array<int, 5>, double> pquv;
      for(int k = 0; k < o1->getnterms(); k++) {
	for(int l = 0; l < o2->getnterms(); l++) {
	  double alpha = o1->getalpha(k), beta = o2->getalpha(l);
	  double zeta = alpha + beta;
	  double Px = (alpha * c1[0] + beta * c2[0]) / zeta,
	    Py = (alpha * c1[1] + beta * c2[1]) / zeta,
	    Pz = (alpha * c1[2] + beta * c2[2]) / zeta;
	  double Kp = M_SQRT2 * std::pow(M_PI, 1.25) / zeta *
	    std::exp(-alpha * beta / zeta *
		     ((c1[0] - c2[0]) * (c1[0] - c2[0]) +
		      (c1[1] - c2[1]) * (c1[1] - c2[1]) +
		      (c1[2] - c2[2]) * (c1[2] - c2[2])));
	  double Dp = o1->getcoef(k) * o1->getnorm(k) *
	    o2->getcoef(l) * o2->getnorm(l);

	  double Rx = Qx - Px,
	    Ry = Qy - Py,
	    Rz = Qz - Pz;
	  double theta2 = zeta * eta / (zeta + eta);
	  double T = theta2 * (Rx * Rx + Ry * Ry + Rz * Rz);

	  double omega = Kp * Kq * Dp * Dq /
	    (std::sqrt(zeta + eta));
	  // Generate the integrals.
	  std::map<int, double> m_ints;
	  for(int m : m_inds) {
	    m_ints[m] = omega * std::pow(2 * theta2, m) *
	      boys_square(m, T);
	  }

	  // Generate the [r] and [p+q]_uv indices.
	  std::map<std::array<int, 4>, double> rm;
	  for(const std::array<int, 5> &ind : pquv_inds) {
	    // Generate [r]
	    std::array<int, 4> ind1 = {
	      ind[0], ind[1], ind[2], 0};
	    double res = compute_rm(ind1, *&rm, m_ints, Rx, Ry, Rz);

	    // Generate (p|q]_uv.
	    if(pquv.count(ind) == 0) {
	      pquv[ind] = 0;
	    }
	    pquv.at(ind) += std::pow(2 * beta, ind[3]) /
	      std::pow(2 * zeta, ind[4]) * res;
	  }
	}
      }

      // Generate (e0|q)_uv
      std::map<std::array<int, 11>, double> apquv;
      for(const std::array<int, 8> &ind : e0quv_inds) {
	std::array<int, 11> ind1 = {
	  ind[0], ind[1], ind[2],
	  0, 0, 0,
	  ind[3], ind[4], ind[5],
	  0, 0};
	double res = compute_apquv(ind1, *&apquv, pquv, c1, c2);

	if(e0quv.count(ind) == 0) {
	  e0quv[ind] = 0;
	}
	e0quv.at(ind) += std::pow(2 * delta, ind[6]) /
	  std::pow(2 * eta, ind[7]) * res;
      }
    }
  }

  // Generate (e0|f0).
  std::map<std::array<int, 11>, double> e0cquv;
  std::map<std::array<int, 6>, double> e0f0;
  for(const std::array<int, 6> &ind : e0f0_inds) {
    std::array<int, 11> ind1 = {
      ind[0], ind[1], ind[2],
      ind[3], ind[4], ind[5],
      0, 0, 0, 0, 0};
    double res = compute_e0cquv(ind1, *&e0cquv, e0quv, c3, c4);
    e0f0[ind] = res;
  }

  // Generate the final integral.
  std::map<std::array<int, 12>, double> abcd;
  return compute_abcd(index, *&abcd, e0f0, c1, c2, c3, c4);
}


double AnalyticIntegral::repulsion(const GaussianOrbital *o1,
				 const GaussianOrbital *o2,
				 const GaussianOrbital *o3,
				 const GaussianOrbital *o4,
				 std::array<double, 3> c1,
				 std::array<double, 3> c2,
				 std::array<double, 3> c3,
				 std::array<double, 3> c4) const {
  double sum = 0;

  for(int i = 0; i < o1->getharms().getsize(); i++) {
    for(int j = 0; j < o2->getharms().getsize(); j++) {
      for(int k = 0; k < o3->getharms().getsize(); k++) {
	for(int l = 0; l < o4->getharms().getsize(); l++) {
	  sum += o1->getharms().getcoef(i) *
	    o2->getharms().getcoef(j) *
	    o3->getharms().getcoef(k) *
	    o4->getharms().getcoef(l) *
	    rep_integral(o1->getharms().gettermorder(i),
			 o2->getharms().gettermorder(j),
			 o3->getharms().gettermorder(k),
			 o4->getharms().gettermorder(l),
			 c1, c2, c3, c4,
			 o1, o2, o3, o4);
	}
      }
    }
  }
  return sum;
}
