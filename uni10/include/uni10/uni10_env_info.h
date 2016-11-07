#ifndef __UNI10_ENV_H__
#define __UNI10_ENV_H__ 

#include <stdio.h>

#include "uni10/uni10_error.h"
#include "uni10/uni10_sys_info.h"

namespace uni10{

  enum exu_type{
    _no    =   0,
    _cpu   =   1,
    _gpu   =   2,
    _mpi   =   3,
    _mgpu  =   4
  };

  class uni10_env{

    public:

      uni10_env(): communicate(false), etype(0){};

      void clear();

      void set_etype(exu_type _etype);

      void set_communicate(uni10_exu_type _communicate);

      void set_memory_info();

      void use_memsize(const uni10_uint64& memsize);

      friend bool uni10_func(load_uni10_rc, _type)(uni10_env& );

      friend void uni10_print_env_info();

    private:

      uni10_exu_type             communicate;
      uni10_exu_type             etype;
      info_type(sysinfo, _type)  uni10_sys_info;

  };

  static uni10_env env_variables;

};

#endif
