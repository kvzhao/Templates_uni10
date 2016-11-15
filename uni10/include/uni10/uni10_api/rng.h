#ifndef __UNI10_RNG_H__
#define __UNI10_RNG_H__

#include <chrono>
#include "uni10/uni10_elem_rng.h"

#define uni10_rng(M, eng, dis, up, dn, seed)\
  do{ \
    uni10_int32  seed1 = seed;\
    uni10_uint64 num = M.elemNum();\
    if(seed1 == -1)\
      seed1 = std::chrono::system_clock::now().time_since_epoch().count();\
    uni10_elem_rng(M, num, eng, dis, up, dn, seed1);\
  }while(0);

#endif
