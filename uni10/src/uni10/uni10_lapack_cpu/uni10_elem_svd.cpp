#include "uni10/uni10_lapack_cpu/uni10_elem_linalg_lapack_cpu.h"

namespace uni10{

  void matrixSVD(const uni10_elem_double64* Mij_ori, uni10_const_bool* isMdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_double64* U, uni10_elem_double64* S, uni10_elem_double64* vT){

    if(!*isMdiag)
      uni10_linalg::matrixSVD(Mij_ori->__elem, *M, *N, U->__elem, S->__elem, vT->__elem);
    else
      uni10_error_msg(true, "%s", "Developping!!!");

  }

  void matrixSVD(const uni10_elem_complex128* Mij_ori, uni10_const_bool* isMdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_complex128* U, uni10_elem_complex128* S, uni10_elem_complex128* vT){

    if(!*isMdiag)
      uni10_linalg::matrixSVD(Mij_ori->__elem, *M, *N, U->__elem, S->__elem, vT->__elem);
    else
      uni10_error_msg(true, "%s", "Developping!!!");

  }

}
