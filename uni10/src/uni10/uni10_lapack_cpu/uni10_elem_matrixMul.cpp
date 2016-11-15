#include "uni10/uni10_lapack_cpu/uni10_elem_linalg_lapack_cpu.h"

namespace uni10{

  void matrixMul(const uni10_elem_double64* A, uni10_const_bool* isAdiag, const uni10_elem_double64* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_double64* C){

    if( !*isAdiag && !*isBdiag ){

      uni10_elem_copy_cpu(C->__elem, A->__elem, A->__elemNum * sizeof(uni10_double64) );
      uni10_linalg::vectorMul(C->__elem, B->__elem, C->__elemNum);

    }
    else if( *isAdiag && !*isBdiag ){

      uni10_uint64 min = std::min(*M, *N);
      for(int i = 0; i < (int)min; i++){
        C->__elem[i * (*N) + i] = A->__elem[i] * B->__elem[i * (*N) + i];
      }
    }
    else if( !*isAdiag && *isBdiag ){

      uni10_uint64 min = std::min(*M, *N);
      for(int i = 0; i < (int)min; i++)
        C->__elem[i * (*N) + i] = B->__elem[i] * A->__elem[i * (*N) + i];

    }
    else{

      uni10_elem_copy_cpu(C->__elem, A->__elem,  C->__elemNum*sizeof(uni10_double64));
      uni10_linalg::vectorMul(C->__elem, B->__elem, C->__elemNum);

    }

  }

  void matrixMul(const uni10_elem_complex128* A, uni10_const_bool* isAdiag, const uni10_elem_complex128* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* C){

    if( !*isAdiag && !*isBdiag ){

      uni10_elem_copy_cpu(C->__elem, A->__elem, A->__elemNum * sizeof(uni10_complex128));
      uni10_linalg::vectorMul(C->__elem, B->__elem, C->__elemNum);

    }
    else if( *isAdiag && !*isBdiag ){

      uni10_uint64 min = std::min(*M, *N);
      for(int i = 0; i < (int)min; i++){
        C->__elem[i * (*N) + i] = A->__elem[i] * B->__elem[i * (*N) + i];
      }

    }
    else if( !*isAdiag && *isBdiag ){

      uni10_uint64 min = std::min(*M, *N);
      for(int i = 0; i < (int)min; i++)
        C->__elem[i * (*N) + i] = B->__elem[i] * A->__elem[i * (*N) + i];

    }
    else{

      uni10_elem_copy_cpu(C->__elem, A->__elem,  A->__elemNum*sizeof(uni10_complex128));
      uni10_linalg::vectorMul(C->__elem, B->__elem, C->__elemNum);

    }

  }

}
