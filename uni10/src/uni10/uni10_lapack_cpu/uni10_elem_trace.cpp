#include "uni10/uni10_lapack_cpu/uni10_elem_linalg_lapack_cpu.h"

namespace uni10{

  void matrixTrace(const uni10_elem_double64* Mij_ori, uni10_const_bool* isMdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_double64* Q, uni10_elem_double64* D, uni10_elem_double64* R){
    
    if(!*isMdiag)
      uni10_linalg::matrixQDR(Mij_ori->__elem, *M, *N, Q->__elem, D->__elem, R->__elem);
    else
      uni10_error_msg(true, "%s", "Developping!!!");

  }

  uni10_complex128 matrixTrace(const uni10_elem_complex128* Mij_ori, uni10_const_bool* isMdiag, const uni10_uint64* M, const uni10_uint64* N){

    uni10_int32 inc = 1;
    uni10_complex128 res = 0.;
    if(!*isMdiag){
      for(uni10_uint64 i = 0; i < min; i++){
        Mij 
      }
    }
    else
      return uni10_linalg::vectorSum(Mij_ori->__elem, *min, inc);

  }

}
