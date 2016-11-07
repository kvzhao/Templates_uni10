// Have to write something.
#include "uni10/uni10_error.h"
#include "uni10/uni10_auxiliary.h"

namespace uni10{

  void uni10_create(){
    
    // Loading the environment setting from .uni10rc or by default.
    if(! load_uni10_rc_cpu(env_variables)){

      // Uni10 is running with CPU only by default.
      env_variables.set_etype(_cpu);
      env_variables.set_communicate(false);

      // Compute the size of the free memory.
      env_variables.set_memory_info();
      //uni10_error_msg(true, "Developping !!!");
    }

  }

  void uni10_destroy(){

    env_variables.clear();
    
  }

  void uni10_print_env_info(){

    env_variables.uni10_sys_info.print_env_info();

  }

#ifdef CPU
  bool load_uni10_rc_cpu(uni10_env& _env){

    bool exsist_rc = true;
    FILE* rcfp = fopen("~/.uni10rc", "r");

    if(!rcfp)
      exsist_rc = false;
    else{

      _env.etype = _cpu;
      _env.communicate = false;

      uni10_error_msg(true, "Developping !!!");

    }

    return exsist_rc;

  }
#endif
#ifdef GPU
  bool load_uni10_rc_gpu(uni10_env& _env){

    //Developping

    bool exsist_rc = true;
    FILE* rcfp = fopen("~/.uni10rc", "r");

    uni10_error_msg(true, "Developping !!!");

    return exsist_rc;

  }
#endif
#ifdef MPI
  bool load_uni10_rc_mpi(uni10_env& _env){

    //Developping

    bool exsist_rc = true;
    FILE* rcfp = fopen("~/.uni10rc", "r");

    uni10_error_msg(true, "Developping !!!");

    return exsist_rc;

  }
#endif

}
