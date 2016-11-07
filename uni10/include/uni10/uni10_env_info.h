#ifndef __UNI10_ENV_H__
#define __UNI10_ENV_H__ 

#include <stdio.h>

#define CPU 1

#include "uni10/uni10_type.h"

namespace uni10{

  enum exu_type{
    exu_no    =   0,
    exu_cpu   =   1,
    exu_gpu   =   2,
    exu_mpi   =   3
  };

  void uni10_create();
  void uni10_destroy();

  class uni10_env{

    public:

      uni10_env(): communicate(false), etype(0), total_memsize(0), free_memsize(0){};

      void clear();

      void set_etype(exu_type _etype);

      void set_communicate(uni10_exu_type _communicate);

      void set_memory_info();

      void use_memsize(const uni10_uint64& memsize);

#ifdef CPU
      friend bool load_uni10_rc_cpu(uni10_env& );
#endif
#ifdef GPU
      friend bool load_uni10_rc_gpu(uni10_env& );
#endif
      friend void uni10_print_env_info();

    private:

      uni10_exu_type    communicate;
      uni10_exu_type    etype;
      uni10_uint64      total_memsize;
      uni10_uint64      free_memsize;

  };


#ifdef CPU
  bool load_uni10_rc_cpu(uni10_env& );
#endif

#ifdef GPU
  bool load_uni10_rc_gpu(uni10_env& );
#endif
  void uni10_print_env_info();

  static uni10_env env_variables;

};

#endif
