#ifndef __UNI10_SYS_INFO_H__
#define __UNI10_SYS_INFO_H__

#include <sys/types.h>
#include <sys/sysctl.h>
#include <sys/sysinfo.h>
#include <unistd.h>

#include "uni10_type.h"
#include "uni10/uni10_error.h"

namespace uni10{

  struct sysinfo_cpu{

    sysinfo_cpu(): total_memsize(0), free_memsize(0), swap_memsize(0), mem_unit(0), core_num(0), status(false){}

    void clear(){
      total_memsize  = 0;
      free_memsize   = 0;
      swap_memsize = 0;
      mem_unit       = 0;
      core_num       = 0;
      status         = false;
    };

    void set_info(){

      struct sysinfo tmp_info;
      sysinfo( & tmp_info );
      total_memsize    = (uni10_uint64)tmp_info.totalram;
      free_memsize     = (uni10_uint64)tmp_info.freeram;
      swap_memsize   = (uni10_uint64)tmp_info.freeswap;
      mem_unit         = (uni10_uint64)tmp_info.mem_unit;
      core_num         = sysconf(_SC_NPROCESSORS_ONLN);
      status = true;

    }

    void print_env_info() const{

      uni10_double64 _total_memsize   = total_memsize  / 1024. / 1024. * mem_unit;
      uni10_double64 _free_memsize    = free_memsize   / 1024. / 1024. * mem_unit;
      uni10_double64 _swap_memsize    = swap_memsize / 1024. / 1024. * mem_unit;
      uni10_uint64 max_len            = floor(log(_total_memsize));

      fprintf(stdout, "\n#######  Uni10 environment information  #######\n");
      fprintf(stdout, "# CPU   cores  : %*.ld \n"  , (int)max_len, core_num);
      fprintf(stdout, "# Total memory : %*.2f MB\n", (int)max_len, _total_memsize);
      fprintf(stdout, "# Free  memory : %*.2f MB\n", (int)max_len, _free_memsize);
      fprintf(stdout, "# Swap  memory : %*.2f MB\n", (int)max_len, _swap_memsize);
      fprintf(stdout, "###############################################\n");
      //fprintf(stdout, "# Memory  unit     : %ld MB\n", env_variables.uni10_sys_info.mem_unit);

    }

    // Parameters
    uni10_uint64 total_memsize;
    uni10_uint64 free_memsize;
    uni10_uint64 swap_memsize;
    uni10_uint64 mem_unit;
    uni10_uint64 core_num;
    bool status;

  };

  struct sysinfo_gpu{
    
    // Developping 
    sysinfo_gpu(){};

    void clear(){};

    void set_info(){};

    void print_env_info(){};

    // Parameters

  };

  struct sysinfo_mpi{

    // Developping 
    sysinfo_mpi();

    void clear(){};

    void set_info(){};

    void print_env_info(){};

    // Parameters
    
  };

}

#endif
