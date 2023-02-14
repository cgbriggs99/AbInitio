/*
 * electron_configs.hpp
 *
 *  Created on: Feb 7, 2023
 *      Author: connor
 */

#ifndef BASIS_ELECTRON_CONFIGS_HPP_
#define BASIS_ELECTRON_CONFIGS_HPP_

#include <initializer_list>
#include <vector>

namespace compchem {

  class GSConfig {
  private :
    int *confs;
    int size_confs;
    int Z;
  public :
    GSConfig(int z, std::initializer_list<int> conf);

    int getShells() const;
    int getZ() const;
    
    int operator[](int index) const;
  };

  const GSConfig &getconfig(int Z);

}



#endif /* BASIS_ELECTRON_CONFIGS_HPP_ */
