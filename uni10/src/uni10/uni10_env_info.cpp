#include "uni10/uni10_error.h"
#include "uni10/uni10_env_info.h"

namespace uni10{

  void uni10_create(){
    
    // Loading the environment setting from .uni10rc or by default.
    if(! load_uni10_rc_cpu(env_variables)){

      //Uni10 is running with CPU only by default.
      env_variables.set_etype(exu_cpu);
      env_variables.set_communicate(false);
      env_variables.set_free_memsize();
      //uni10_error_msg(true, "Developping !!!");

    }

    // Compute the size of the free memory.
    
  }

  void uni10_destroy(){

    env_variables.clear();
    
  }

  void uni10_print_env_info(){

      uni10_error_msg(true, "Developping !!!");
    
  }

  bool load_uni10_rc_cpu(uni10_env& _env){

    bool exsist_rc = true;
    FILE* rcfp = fopen("~/.uni10rc", "r");

    if(!rcfp)
      exsist_rc = false;
    else{

      _env.etype = exu_cpu;
      _env.communicate = false;

      uni10_error_msg(true, "Developping !!!");

    }
    
    return exsist_rc;

  }

  void uni10_env::set_free_memsize(){

      uni10_error_msg(true, "Developping !!!");

  }

}
