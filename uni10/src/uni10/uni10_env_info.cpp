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

    uni10_sys_info.set_info();

  }

  void uni10_env::clear(){

    communicate   = false;
    etype         = _no;
    uni10_sys_info.clear();

  }

  void uni10_env::use_memsize(const uni10_uint64& memsize){
    env_variables.uni10_sys_info.free_memsize -= memsize;
  }

}
