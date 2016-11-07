#include "uni10/uni10_error.h"
#include "uni10/uni10_env_info.h"

namespace uni10{

  void uni10_env::set_etype(exu_type _etype){
    etype = _etype;
  }

  void uni10_env::set_communicate(uni10_exu_type _communicate){
    communicate = _communicate;
  };

  void uni10_env::set_memory_info(){
      uni10_error_msg(true, "Developping !!!");
  }

  void uni10_env::use_memsize(const uni10_uint64& memsize){
      free_memsize -= memsize;
  }

  void uni10_env::clear(){

    communicate   = false;
    etype         = exu_no;
    total_memsize = 0;
    free_memsize  = 0;
    
  }

  void uni10_create(){
    
    // Loading the environment setting from .uni10rc or by default.
    if(! load_uni10_rc_cpu(env_variables)){

      //Uni10 is running with CPU only by default.
      env_variables.set_etype(exu_cpu);
      env_variables.set_communicate(false);
      env_variables.set_memory_info();
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

}
