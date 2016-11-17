#include "uni10/uni10_lapack_cpu/uni10_elem_linalg_lapack_cpu.h"

namespace uni10{

  void matrixSub(uni10_elem_double64* A, uni10_bool* isAdiag, const uni10_elem_double64* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, const uni10_uint64* N){

    if(isAdiag && !isBdiag){

      uni10_uint64 elemNum = A->__elemNum;

      uni10_double64* _elem = (uni10_double64*)malloc(A->__elemNum * sizeof(uni10_double64));

      uni10_elem_copy_cpu(_elem, A->__elem, A->__elemNum*sizeof(uni10_double64));

      *isAdiag = false;

      A->init(*M, *N, false, B->__elem);

      uni10_double64 alpha = -1.;

      uni10_linalg::vectorScal(alpha , A->__elem, A->__elemNum);

      for(int i = 0; i < (int)elemNum; i++)
        A->__elem[i*(*N)+i] += _elem[i];

      free(_elem);

    }
    else if(!isAdiag && isBdiag){

      uni10_uint64 elemNum = B->__elemNum;

      for(int i = 0; i < (int)elemNum; i++)
        A->__elem[i*(*N)+i] -= B->__elem[i];

    }
    else
      uni10_linalg::vectorSub(A->__elem, B->__elem, B->__elemNum);


  }

  void matrixSub(uni10_elem_complex128* A, uni10_bool* isAdiag, const uni10_elem_complex128* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, const uni10_uint64* N){

    if(isAdiag && !isBdiag){

      uni10_uint64 elemNum = A->__elemNum;

      uni10_complex128* _elem = (uni10_complex128*)malloc(A->__elemNum * sizeof(uni10_complex128));

      uni10_elem_copy_cpu(_elem, A->__elem, A->__elemNum*sizeof(uni10_complex128));

      *isAdiag = false;

      A->init(*M, *N, false, B->__elem);

      uni10_double64 alpha = -1.;

      uni10_linalg::vectorScal(alpha , A->__elem, A->__elemNum);

      for(int i = 0; i < (int)elemNum; i++)
        A->__elem[i*(*N)+i] += _elem[i];

      free(_elem);

    }
    else if(!isAdiag && isBdiag){

      uni10_uint64 elemNum = B->__elemNum;

      for(int i = 0; i < (int)elemNum; i++)
        A->__elem[i*(*N)+i] -= B->__elem[i];

    }
    else
      uni10_linalg::vectorSub(A->__elem, B->__elem, B->__elemNum);

  }

  void matrixSub(uni10_elem_complex128* A, uni10_bool* isAdiag, const uni10_elem_double64* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, const uni10_uint64* N){

    if(isAdiag && !isBdiag){

      uni10_uint64 elemNum = A->__elemNum;

      uni10_complex128* _elem = (uni10_complex128*)malloc(A->__elemNum * sizeof(uni10_complex128));

      uni10_elem_copy_cpu(_elem, A->__elem, A->__elemNum*sizeof(uni10_double64));

      *isAdiag = false;

      A->init(*M, *N, false, NULL);

      uni10_elem_cast_cpu(A->__elem, B->__elem, B->__elemNum);

      uni10_double64 alpha = -1.;

      uni10_linalg::vectorScal(alpha , A->__elem, A->__elemNum);

      for(int i = 0; i < (int)elemNum; i++)
        A->__elem[i*(*N)+i] += _elem[i];

      free(_elem);

    }
    else if(!isAdiag && isBdiag){

      uni10_uint64 elemNum = B->__elemNum;

      for(int i = 0; i < (int)elemNum; i++)
        A->__elem[i*(*N)+i] -= B->__elem[i];

    }
    else
      uni10_linalg::vectorSub(A->__elem, B->__elem, B->__elemNum);

  }

  void matrixSub(const uni10_elem_double64* A, uni10_const_bool* isAdiag, const uni10_elem_double64* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_double64* C){

    if( !*isAdiag && !*isBdiag ){

      uni10_elem_copy_cpu(C->__elem, A->__elem, A->__elemNum * sizeof(uni10_double64) );
      uni10_linalg::vectorSub(C->__elem, B->__elem, C->__elemNum);

    }
    else if( *isAdiag && !*isBdiag ){

      uni10_elem_copy_cpu(C->__elem, B->__elem, B->__elemNum * sizeof(uni10_double64) );
      uni10_uint64 min = std::min(*M, *N);
      uni10_linalg::vectorScal(-1., C->__elem, C->__elemNum);
      for(int i = 0; i < (int)min; i++){
        C->__elem[i * (*N) + i] += A->__elem[i];
      }
    }
    else if( !*isAdiag && *isBdiag ){

      uni10_elem_copy_cpu(C->__elem, A->__elem, A->__elemNum * sizeof(uni10_double64) );
      uni10_uint64 min = std::min(*M, *N);
      for(int i = 0; i < (int)min; i++)
        C->__elem[i * (*N) + i] -= B->__elem[i];

    }
    else{

      uni10_elem_copy_cpu(C->__elem, A->__elem,  C->__elemNum*sizeof(uni10_double64));
      uni10_linalg::vectorSub(C->__elem, B->__elem, C->__elemNum);

    }

  }

  void matrixSub(const uni10_elem_complex128* A, uni10_const_bool* isAdiag, const uni10_elem_complex128* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* C){

    if( !*isAdiag && !*isBdiag ){

      uni10_elem_copy_cpu(C->__elem, A->__elem, A->__elemNum * sizeof(uni10_complex128));
      uni10_linalg::vectorSub(C->__elem, B->__elem, C->__elemNum);

    }
    else if( *isAdiag && !*isBdiag ){

      uni10_elem_copy_cpu(C->__elem, B->__elem, B->__elemNum * sizeof(uni10_complex128));
      uni10_uint64 min = std::min(*M, *N);
      uni10_linalg::vectorScal(-1., C->__elem, C->__elemNum);
      for(int i = 0; i < (int)min; i++){
        C->__elem[i * (*N) + i] += A->__elem[i];
      }

    }
    else if( !*isAdiag && *isBdiag ){

      uni10_elem_copy_cpu(C->__elem, A->__elem, A->__elemNum * sizeof(uni10_complex128));
      uni10_uint64 min = std::min(*M, *N);
      for(int i = 0; i < (int)min; i++)
        C->__elem[i * (*N) + i] -= B->__elem[i];

    }
    else{

      uni10_elem_copy_cpu(C->__elem, A->__elem,  A->__elemNum*sizeof(uni10_complex128));
      uni10_linalg::vectorSub(C->__elem, B->__elem, C->__elemNum);

    }

  }

  void matrixSub(const uni10_elem_double64* A, uni10_const_bool* isAdiag, const uni10_elem_complex128* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* C){

    uni10_error_msg(true, "%s", "Debugging !!!");

    if( !*isAdiag && !*isBdiag ){

      uni10_elem_copy_cpu(C->__elem, B->__elem, B->__elemNum * sizeof(uni10_complex128) );
      uni10_linalg::vectorSub(C->__elem, A->__elem, C->__elemNum);

    }
    else if( *isAdiag && !*isBdiag ){

      uni10_elem_copy_cpu(C->__elem, B->__elem, B->__elemNum * sizeof(uni10_complex128) );
      uni10_uint64 min = std::min(*M, *N);
      for(int i = 0; i < (int)min; i++)
        C->__elem[i * (*N) + i] -= A->__elem[i];

    }
    else if( !*isAdiag && *isBdiag ){

      uni10_elem_cast_cpu(C->__elem, A->__elem, C->__elemNum);
      uni10_uint64 min = std::min(*M, *N);
      for(int i = 0; i < (int)min; i++)
        C->__elem[i * (*N) + i] -= B->__elem[i];

    }
    else{

      uni10_elem_copy_cpu(C->__elem, B->__elem, C->__elemNum*sizeof(uni10_complex128));
      uni10_linalg::vectorSub(C->__elem, A->__elem, C->__elemNum);

    }

  }

  void matrixSub(const uni10_elem_complex128* A, uni10_const_bool* isAdiag, const uni10_elem_double64* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* C){

    uni10_error_msg(true, "%s", "Debugging !!!");

    if( !*isAdiag && !*isBdiag ){

      uni10_elem_copy_cpu(C->__elem, A->__elem, A->__elemNum * sizeof(uni10_complex128) );
      uni10_linalg::vectorSub(C->__elem, B->__elem, C->__elemNum);

    }
    else if( *isAdiag && !*isBdiag ){

      uni10_elem_cast_cpu(C->__elem, B->__elem, B->__elemNum);
      uni10_uint64 min = std::min(*M, *N);
      for(int i = 0; i < (int)min; i++)
        C->__elem[i * (*N) + i] -= A->__elem[i];

    }
    else if( !*isAdiag && *isBdiag ){

      uni10_elem_copy_cpu(C->__elem, A->__elem, A->__elemNum * sizeof(uni10_complex128) );
      uni10_uint64 min = std::min(*M, *N);
      for(int i = 0; i < (int)min; i++)
        C->__elem[i * (*N) + i] -= B->__elem[i];

    }
    else{

      uni10_elem_copy_cpu(C->__elem, A->__elem, A->__elemNum*sizeof(uni10_complex128));
      uni10_linalg::vectorSub(C->__elem, A->__elem, C->__elemNum);

    }

  }

}
