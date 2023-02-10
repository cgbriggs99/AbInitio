 
#ifndef __COPYABLE_HPP__
#define __COPYABLE_HPP__

namespace compchem {

  class Copyable {

  public :
    virtual ~Copyable() = 0;

    virtual Copyable &copy() = 0;
  };

}

#endif
