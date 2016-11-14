#include "uni10/uni10_lapack_cpu/uni10_tools_cpu.h"
#include "uni10/uni10_lapack_cpu/uni10_elem_linalg_lapack_cpu.h"

namespace uni10{

  // BLAS
  //
  // UNI10_DOUBLE64
  void vectorAdd(uni10_elem_double64* Y, const uni10_elem_double64* X, const uni10_uint64* N){

    uni10_linalg::vectorAdd(Y->__elem, X->__elem, *N);

  }

  void vectorSub(uni10_elem_double64* Y, const uni10_elem_double64* X, const uni10_uint64* N){

    uni10_linalg::vectorSub(Y->__elem, X->__elem, *N);

  }

  void vectorMul(uni10_elem_double64* Y, const uni10_elem_double64* X, const uni10_uint64* N){

    uni10_linalg::vectorSub(Y->__elem, X->__elem, *N);

  }

  void vectorScal(uni10_double64* a, uni10_elem_double64* X, uni10_uint64* N){

    uni10_linalg::vectorScal(*a, X->__elem, *N);

  }

  void vectorExp(uni10_double64* a, uni10_elem_double64* X, uni10_uint64* N){

    uni10_linalg::vectorExp(*a, X->__elem, *N);

  }

  uni10_double64 vectorSum(const uni10_elem_double64* X, const uni10_uint64* N, uni10_int32* inc){

    return uni10_linalg::vectorSum(X->__elem, *N, *inc);

  }

  uni10_double64 vectorNorm(const uni10_elem_double64* X, const uni10_uint64* N, uni10_int32* inc){

    return uni10_linalg::vectorNorm(X->__elem, *N, *inc);

  }
  
  void matrixAdd(const uni10_elem_double64* A, uni10_const_bool* isAdiag, const uni10_elem_double64* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_double64* C){

    if( !*isAdiag && !*isBdiag ){

      uni10_elem_copy_cpu(C->__elem, B->__elem, B->__elemNum * sizeof(uni10_double64) );
      uni10_linalg::vectorAdd(C->__elem, A->__elem, C->__elemNum);

    }
    else if( *isAdiag && !*isBdiag ){

      uni10_elem_copy_cpu(C->__elem, B->__elem, B->__elemNum * sizeof(uni10_double64) );
      uni10_uint64 min = min(*M, *N);
      for(int i = 0; i < (int)min; i++)
        C->__elem[i * (*N) + i] += A->__elem[i];

    }
    else if( !*isAdiag && *isBdiag ){

      uni10_elem_copy_cpu(C->__elem, A->__elem, A->__elemNum * sizeof(uni10_double64) );
      uni10_uint64 min = min(*M, *N);
      for(int i = 0; i < (int)min; i++)
        C->__elem[i * (*N) + i] += B->__elem[i];

    }
    else{

      uni10_elem_copy_cpu(C->__elem, B->__elem,  C->__elemNum*sizeof(uni10_double64));
      uni10_linalg::vectorAdd(C->__elem, A->__elem, C->__elemNum);

    }

  }

  void matrixSub(const uni10_elem_double64* A, uni10_const_bool* isAdiag, const uni10_elem_double64* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_double64* C){

    if( !*isAdiag && !*isBdiag ){

      uni10_elem_copy_cpu(C->__elem, A->__elem, A->__elemNum * sizeof(uni10_double64) );
      uni10_linalg::vectorSub(C->__elem, B->__elem, C->__elemNum);

    }
    else if( *isAdiag && !*isBdiag ){

      uni10_elem_copy_cpu(C->__elem, B->__elem, B->__elemNum * sizeof(uni10_double64) );
      uni10_uint64 min = min(*M, *N);
      uni10_linalg::vectorScal(-1., C->__elem, C->__elemNum);
      for(int i = 0; i < (int)min; i++){
        C->__elem[i * (*N) + i] += A->__elem[i];
      }
    }
    else if( !*isAdiag && *isBdiag ){

      uni10_elem_copy_cpu(C->__elem, A->__elem, A->__elemNum * sizeof(uni10_double64) );
      uni10_uint64 min = min(*M, *N);
      for(int i = 0; i < (int)min; i++)
        C->__elem[i * (*N) + i] -= B->__elem[i];

    }
    else{

      uni10_elem_copy_cpu(C->__elem, A->__elem,  C->__elemNum*sizeof(uni10_double64));
      uni10_linalg::vectorSub(C->__elem, B->__elem, C->__elemNum);

    }

  }

  void matrixMul(const uni10_elem_double64* A, uni10_const_bool* isAdiag, const uni10_elem_double64* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_double64* C){

    if( !*isAdiag && !*isBdiag ){

      uni10_elem_copy_cpu(C->__elem, A->__elem, A->__elemNum * sizeof(uni10_double64) );
      uni10_linalg::vectorMul(C->__elem, B->__elem, C->__elemNum);

    }
    else if( *isAdiag && !*isBdiag ){

      uni10_uint64 min = min(*M, *N);
      for(int i = 0; i < (int)min; i++){
        C->__elem[i * (*N) + i] = A->__elem[i] * B->__elem[i * (*N) + i];
      }
    }
    else if( !*isAdiag && *isBdiag ){

      uni10_uint64 min = min(*M, *N);
      for(int i = 0; i < (int)min; i++)
        C->__elem[i * (*N) + i] = B->__elem[i] * A->__elem[i * (*N) + i];

    }
    else{

      uni10_elem_copy_cpu(C->__elem, A->__elem,  C->__elemNum*sizeof(uni10_double64));
      uni10_linalg::vectorMul(C->__elem, B->__elem, C->__elemNum);

    }

  }

  void matrixDot(const uni10_elem_double64* A, uni10_const_bool* isAdiag, const uni10_elem_double64* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, const uni10_uint64* N, const uni10_uint64* K, uni10_elem_double64* C){

    if( !*isAdiag && !*isBdiag ){

      uni10_linalg::matrixDot(A->__elem, B->__elem, *M, *N, *K, C->__elem);

    }
    else if( *isAdiag && !*isBdiag ){

      uni10_elem_copy_cpu(C->__elem, B->__elem, B->__elemNum*sizeof(uni10_double64));
      uni10_linalg::diagRowMul(C->__elem, A->__elem, min(*M, *K), *N);

    }
    else if( !*isAdiag && *isBdiag ){

      uni10_uint64 data_col = min(*K, *N);

      for(int r = 0; r < (int)*M; r++)
        uni10_elem_copy_cpu(&C->__elem[r*(*N)], &A->__elem[r*(*K)], data_col*sizeof(uni10_double64));

      uni10_linalg::diagColMul(C->__elem, A->__elem, *M, data_col);

    }
    else{

      uni10_uint64 min = min(A->__elemNum, B->__elemNum);
      uni10_elem_copy_cpu(C->__elem, A->__elem, min*sizeof(uni10_double64));

      uni10_linalg::vectorMul(C->__elem, B->__elem, min);

    }

  }

  void setTranspose(const uni10_elem_double64* A, const uni10_uint64* M, const uni10_uint64* N, uni10_elem_double64* AT){

    uni10_linalg::setTranspose(A->__elem, *M, *N, AT->__elem);

  }

  void setTranspose(uni10_elem_double64* A, uni10_uint64* M, uni10_uint64* N){
    
    uni10_linalg::setTranspose(A->__elem, *M, *N);

  }

  void setDagger(const uni10_elem_double64* A, const uni10_uint64* M, const uni10_uint64* N, uni10_elem_double64* AT){

    uni10_linalg::setDagger(A->__elem, *M, *N, AT->__elem);

  }

  void setDagger(uni10_elem_double64* A, uni10_uint64* M, uni10_uint64* N){

    uni10_linalg::setDagger(A->__elem, *M, *N);

  }

  void setConjugate(const uni10_elem_double64* A, const uni10_uint64* N, uni10_elem_double64* A_conj){

    uni10_elem_copy_cpu(A_conj->__elem, A->__elem, *N *sizeof(uni10_elem_double64));

  }

  void setConjugate(uni10_elem_double64* A, uni10_uint64* N){

    //For function overload. Nothing to do.
    //
  }

  // BLAS
  //
  // UNI10_COMPLEX64
  void vectorAdd(uni10_elem_complex128* Y, const uni10_elem_complex128* X, const uni10_uint64* N){

    uni10_linalg::vectorAdd(Y->__elem, X->__elem, *N);

  }

  void vectorSub(uni10_elem_complex128* Y, const uni10_elem_complex128* X, const uni10_uint64* N){

    uni10_linalg::vectorSub(Y->__elem, X->__elem, *N);

  }

  void vectorMul(uni10_elem_complex128* Y, const uni10_elem_complex128* X, const uni10_uint64* N){

    uni10_linalg::vectorMul(Y->__elem, X->__elem, *N);

  }  

  void vectorScal(uni10_complex128* a, uni10_elem_complex128* X, uni10_uint64* N){

    uni10_linalg::vectorScal(*a, X->__elem, *N);

  }

  void vectorExp(uni10_complex128* a, uni10_elem_complex128* X, uni10_uint64* N){

    uni10_linalg::vectorExp(*a, X->__elem, *N);

  }
  
  void matrixAdd(const uni10_elem_complex128* A, uni10_const_bool* isAdiag, const uni10_elem_complex128* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* C){

    if( !*isAdiag && !*isBdiag ){

      uni10_elem_copy_cpu(C->__elem, B->__elem, B->__elemNum * sizeof(uni10_complex128));
      uni10_linalg::vectorAdd(C->__elem, A->__elem, C->__elemNum);

    }
    else if( *isAdiag && !*isBdiag ){

      uni10_elem_copy_cpu(C->__elem, B->__elem, B->__elemNum * sizeof(uni10_complex128));
      uni10_uint64 min = min(*M, *N);
      for(int i = 0; i < (int)min; i++)
        C->__elem[i * (*N) + i] += A->__elem[i];

    }
    else if( !*isAdiag && *isBdiag ){

      uni10_elem_copy_cpu(C->__elem, A->__elem, A->__elemNum * sizeof(uni10_complex128));
      uni10_uint64 min = min(*M, *N);
      for(int i = 0; i < (int)min; i++)
        C->__elem[i * (*N) + i] += B->__elem[i];

    }
    else{

      uni10_elem_copy_cpu(C->__elem, B->__elem,  C->__elemNum*sizeof(uni10_complex128));
      uni10_linalg::vectorAdd(C->__elem, A->__elem, C->__elemNum);

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
      uni10_uint64 min = min(*M, *N);
      uni10_linalg::vectorScal(-1., C->__elem, C->__elemNum);
      for(int i = 0; i < (int)min; i++){
        C->__elem[i * (*N) + i] += A->__elem[i];
      }

    }
    else if( !*isAdiag && *isBdiag ){

      uni10_elem_copy_cpu(C->__elem, A->__elem, A->__elemNum * sizeof(uni10_complex128));
      uni10_uint64 min = min(*M, *N);
      for(int i = 0; i < (int)min; i++)
        C->__elem[i * (*N) + i] -= B->__elem[i];

    }
    else{

      uni10_elem_copy_cpu(C->__elem, A->__elem,  A->__elemNum*sizeof(uni10_complex128));
      uni10_linalg::vectorSub(C->__elem, B->__elem, C->__elemNum);

    }

  }

  void matrixMul(const uni10_elem_complex128* A, uni10_const_bool* isAdiag, const uni10_elem_complex128* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* C){

    if( !*isAdiag && !*isBdiag ){

      uni10_elem_copy_cpu(C->__elem, A->__elem, A->__elemNum * sizeof(uni10_complex128));
      uni10_linalg::vectorMul(C->__elem, B->__elem, C->__elemNum);

    }
    else if( *isAdiag && !*isBdiag ){

      uni10_uint64 min = min(*M, *N);
      for(int i = 0; i < (int)min; i++){
        C->__elem[i * (*N) + i] = A->__elem[i] * B->__elem[i * (*N) + i];
      }

    }
    else if( !*isAdiag && *isBdiag ){

      uni10_uint64 min = min(*M, *N);
      for(int i = 0; i < (int)min; i++)
        C->__elem[i * (*N) + i] = B->__elem[i] * A->__elem[i * (*N) + i];

    }
    else{

      uni10_elem_copy_cpu(C->__elem, A->__elem,  A->__elemNum*sizeof(uni10_complex128));
      uni10_linalg::vectorMul(C->__elem, B->__elem, C->__elemNum);

    }

  }

  void matrixDot(const uni10_elem_complex128* A, uni10_const_bool* isAdiag, const uni10_elem_complex128* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, const uni10_uint64* N, const uni10_uint64* K, uni10_elem_complex128* C){

    if( !*isAdiag && !*isBdiag ){

      uni10_linalg::matrixDot(A->__elem, B->__elem, *M, *N, *K, C->__elem);

    }
    else if( *isAdiag && !*isBdiag ){

      uni10_elem_copy_cpu(C->__elem, B->__elem, B->__elemNum*sizeof(uni10_complex128));
      uni10_linalg::diagRowMul(C->__elem, A->__elem, min(*M, *K), *N);

    }
    else if( !*isAdiag && *isBdiag ){

      uni10_uint64 data_col = min(*K, *N);

      for(int r = 0; r < (int)*M; r++)
        uni10_elem_copy_cpu(&C->__elem[r*(*N)], &A->__elem[r*(*K)], data_col*sizeof(uni10_complex128));

      uni10_linalg::diagColMul(C->__elem, A->__elem, *M, data_col);

    }
    else{

      uni10_uint64 min = min(A->__elemNum, B->__elemNum);
      uni10_elem_copy_cpu(C->__elem, A->__elem, min*sizeof(uni10_complex128));

      uni10_linalg::vectorMul(C->__elem, B->__elem, min);

    }

  }

  void setTranspose(const uni10_elem_complex128* A, const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* AT){

    uni10_linalg::setTranspose(A->__elem, *M, *N, AT->__elem);

  }

  void setTranspose(uni10_elem_complex128* A, uni10_uint64* M, uni10_uint64* N){

    uni10_linalg::setTranspose(A->__elem, *M, *N);

  }

  void setDagger(const uni10_elem_complex128* A, const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* AT){

    uni10_linalg::setDagger(A->__elem, *M, *N, AT->__elem);

  }

  void setDagger(uni10_elem_complex128* A, uni10_uint64* M, uni10_uint64* N){

    uni10_linalg::setDagger(A->__elem, *M, *N);

  }

  void setConjugate(const uni10_elem_complex128* A, const uni10_uint64* N, uni10_elem_complex128* A_conj){

    uni10_linalg::setConjugate(A->__elem, *N, A_conj->__elem);
    //For function overload. Nothing to do.
    //
  }

  void setConjugate(uni10_elem_complex128* A, uni10_uint64* N){

    uni10_linalg::setConjugate(A->__elem, *N);

  }

  // BLAS
  //
  // MIX
  //
  void vectorAdd(uni10_elem_complex128* Y, const uni10_elem_double64* X, const uni10_uint64* N){

    uni10_linalg::vectorAdd(Y->__elem, X->__elem, *N);

  }

  void vectorSub(uni10_elem_complex128* Y, const uni10_elem_double64* X, const uni10_uint64* N){

    uni10_linalg::vectorSub(Y->__elem, X->__elem, *N);

  }

  void vectorMul(uni10_elem_complex128* Y, const uni10_elem_double64* X, const uni10_uint64* N){

    uni10_linalg::vectorMul(Y->__elem, X->__elem, *N);

  }

  void vectorScal(uni10_double64* a, uni10_elem_complex128* X, uni10_uint64* N){

    uni10_linalg::vectorScal(*a, X->__elem, *N);

  }

  void vectorExp(uni10_double64* a, uni10_elem_complex128* X, uni10_uint64* N){

    uni10_linalg::vectorExp(*a, X->__elem, *N);

  }

  void matrixAdd(const uni10_elem_double64* A, uni10_const_bool* isAdiag, const uni10_elem_complex128* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* C){

    if( !*isAdiag && !*isBdiag ){

      uni10_elem_copy_cpu(C->__elem, B->__elem, B->__elemNum * sizeof(uni10_complex128) );
      uni10_linalg::vectorAdd(C->__elem, A->__elem, C->__elemNum);

    }
    else if( *isAdiag && !*isBdiag ){

      uni10_elem_copy_cpu(C->__elem, B->__elem, B->__elemNum * sizeof(uni10_complex128) );
      uni10_uint64 min = min(*M, *N);
      for(int i = 0; i < (int)min; i++)
        C->__elem[i * (*N) + i] += A->__elem[i];

    }
    else if( !*isAdiag && *isBdiag ){

      uni10_elem_cast_cpu(C->__elem, A->__elem, C->__elemNum);
      uni10_uint64 min = min(*M, *N);
      for(int i = 0; i < (int)min; i++)
        C->__elem[i * (*N) + i] += B->__elem[i];

    }
    else{

      uni10_elem_copy_cpu(C->__elem, B->__elem,  C->__elemNum*sizeof(uni10_complex128));
      uni10_linalg::vectorAdd(C->__elem, A->__elem, C->__elemNum);

    }

  }

  void matrixSub(const uni10_elem_double64* A, uni10_const_bool* isAdiag, const uni10_elem_complex128* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* C){

    uni10_error_msg(true, "Debugging !!!");

    if( !*isAdiag && !*isBdiag ){

      uni10_elem_copy_cpu(C->__elem, B->__elem, B->__elemNum * sizeof(uni10_complex128) );
      uni10_linalg::vectorSub(C->__elem, A->__elem, C->__elemNum);

    }
    else if( *isAdiag && !*isBdiag ){

      uni10_elem_copy_cpu(C->__elem, B->__elem, B->__elemNum * sizeof(uni10_complex128) );
      uni10_uint64 min = min(*M, *N);
      for(int i = 0; i < (int)min; i++)
        C->__elem[i * (*N) + i] -= A->__elem[i];

    }
    else if( !*isAdiag && *isBdiag ){

      uni10_elem_cast_cpu(C->__elem, A->__elem, C->__elemNum);
      uni10_uint64 min = min(*M, *N);
      for(int i = 0; i < (int)min; i++)
        C->__elem[i * (*N) + i] -= B->__elem[i];

    }
    else{

      uni10_elem_copy_cpu(C->__elem, B->__elem, C->__elemNum*sizeof(uni10_complex128));
      uni10_linalg::vectorSub(C->__elem, A->__elem, C->__elemNum);

    }

  }

  void matrixDot(const uni10_elem_double64* A, uni10_const_bool* isAdiag, const uni10_elem_complex128* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, const uni10_uint64* N, const uni10_uint64* K, uni10_elem_complex128* C){

    if( !*isAdiag && !*isBdiag ){

      uni10_linalg::matrixDot(A->__elem, B->__elem, *M, *N, *K, C->__elem);

    }
    else if( *isAdiag && !*isBdiag ){

      uni10_elem_copy_cpu(C->__elem, B->__elem, B->__elemNum*sizeof(uni10_complex128));
      uni10_linalg::diagRowMul(C->__elem, A->__elem, min(*M, *K), *N);

    }
    else if( !*isAdiag && *isBdiag ){

      uni10_uint64 data_col = min(*K, *N);

      for(int r = 0; r < (int)*M; r++)
        uni10_elem_cast_cpu(&C->__elem[r*(*N)], &A->__elem[r*(*K)], data_col);

      uni10_linalg::diagColMul(C->__elem, A->__elem, *M, data_col);

    }
    else{

      uni10_uint64 min = min(A->__elemNum, B->__elemNum);
      uni10_elem_copy_cpu(C->__elem, A->__elem, min*sizeof(uni10_complex128));

      uni10_linalg::vectorMul(C->__elem, B->__elem, min);

    }

  }

  void matrixAdd(const uni10_elem_complex128* A, uni10_const_bool* isAdiag, const uni10_elem_double64* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* C){

    if( !*isAdiag && !*isBdiag ){

      uni10_elem_copy_cpu(C->__elem, A->__elem, A->__elemNum * sizeof(uni10_complex128) );
      uni10_linalg::vectorAdd(C->__elem, B->__elem, C->__elemNum);

    }
    else if( *isAdiag && !*isBdiag ){

      uni10_elem_cast_cpu(C->__elem, B->__elem, B->__elemNum);
      uni10_uint64 min = min(*M, *N);
      for(int i = 0; i < (int)min; i++)
        C->__elem[i * (*N) + i] += A->__elem[i];

    }
    else if( !*isAdiag && *isBdiag ){

      uni10_elem_copy_cpu(C->__elem, A->__elem, A->__elemNum * sizeof(uni10_complex128) );
      uni10_uint64 min = min(*M, *N);
      for(int i = 0; i < (int)min; i++)
        C->__elem[i * (*N) + i] += B->__elem[i];

    }
    else{

      uni10_elem_copy_cpu(C->__elem, A->__elem, A->__elemNum*sizeof(uni10_complex128));
      uni10_linalg::vectorAdd(C->__elem, A->__elem, C->__elemNum);

    }

  }

  
  void matrixSub(const uni10_elem_complex128* A, uni10_const_bool* isAdiag, const uni10_elem_double64* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* C){

    uni10_error_msg(true, "Debugging !!!");

    if( !*isAdiag && !*isBdiag ){

      uni10_elem_copy_cpu(C->__elem, A->__elem, A->__elemNum * sizeof(uni10_complex128) );
      uni10_linalg::vectorSub(C->__elem, B->__elem, C->__elemNum);

    }
    else if( *isAdiag && !*isBdiag ){

      uni10_elem_cast_cpu(C->__elem, B->__elem, B->__elemNum);
      uni10_uint64 min = min(*M, *N);
      for(int i = 0; i < (int)min; i++)
        C->__elem[i * (*N) + i] -= A->__elem[i];

    }
    else if( !*isAdiag && *isBdiag ){

      uni10_elem_copy_cpu(C->__elem, A->__elem, A->__elemNum * sizeof(uni10_complex128) );
      uni10_uint64 min = min(*M, *N);
      for(int i = 0; i < (int)min; i++)
        C->__elem[i * (*N) + i] -= B->__elem[i];

    }
    else{

      uni10_elem_copy_cpu(C->__elem, A->__elem, A->__elemNum*sizeof(uni10_complex128));
      uni10_linalg::vectorSub(C->__elem, A->__elem, C->__elemNum);

    }

  }

  void matrixDot(const uni10_elem_complex128* A, uni10_const_bool* isAdiag, const uni10_elem_double64* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, const uni10_uint64* N, const uni10_uint64* K, uni10_elem_complex128* C){

    if( !*isAdiag && !*isBdiag ){

      uni10_linalg::matrixDot(A->__elem, B->__elem, *M, *N, *K, C->__elem);

    }
    else if( *isAdiag && !*isBdiag ){

      uni10_elem_cast_cpu(C->__elem, B->__elem, B->__elemNum);
      uni10_linalg::diagRowMul(C->__elem, A->__elem, min(*M, *K), *N);

    }
    else if( !*isAdiag && *isBdiag ){

      uni10_uint64 data_col = min(*K, *N);

      for(int r = 0; r < (int)*M; r++)
        uni10_elem_copy_cpu(&C->__elem[r*(*N)], &A->__elem[r*(*K)], data_col*sizeof(uni10_complex128));

      uni10_linalg::diagColMul(C->__elem, A->__elem, *M, data_col);

    }
    else{

      uni10_uint64 min = min(A->__elemNum, B->__elemNum);
      uni10_elem_copy_cpu(C->__elem, A->__elem, min*sizeof(uni10_complex128));

      uni10_linalg::vectorMul(C->__elem, B->__elem, min);

    }

  }

  uni10_complex128 vectorSum(const uni10_elem_complex128* X, const uni10_uint64* N, uni10_int32* inc){

    return uni10_linalg::vectorSum(X->__elem, *N, *inc);
  }


  uni10_double64 vectorNorm(const uni10_elem_complex128* X, const uni10_uint64* N, uni10_int32* inc){

    return uni10_linalg::vectorNorm(X->__elem, *N, *inc);

  }

  // LAPACK
  //
  // UNI10_DOUBLE64
  void matrixQR(const uni10_elem_double64* Mij_ori, uni10_const_bool* isMdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_double64* Q, uni10_elem_double64* R){
    
    if(!*isMdiag)
      uni10_linalg::matrixQR(Mij_ori->__elem, *M, *N, Q->__elem, R->__elem);
    else
      uni10_error_msg(true, "Developping!!!");

  }

  void matrixRQ(const uni10_elem_double64* Mij_ori, uni10_const_bool* isMdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_double64* Q, uni10_elem_double64* R){

    uni10_error_msg(true, "Developping!!!(Have bugs)");
    if(!*isMdiag)
      uni10_linalg::matrixRQ(Mij_ori->__elem, *M, *N, Q->__elem, R->__elem);
    else
      uni10_error_msg(true, "Developping!!!");

  }

  void matrixQL(const uni10_elem_double64* Mij_ori, uni10_const_bool* isMdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_double64* Q, uni10_elem_double64* L){

    if(!*isMdiag)
      uni10_linalg::matrixQL(Mij_ori->__elem, *M, *N, Q->__elem, L->__elem);
    else
      uni10_error_msg(true, "Developping!!!");

  }

  void matrixLQ(const uni10_elem_double64* Mij_ori, uni10_const_bool* isMdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_double64* Q, uni10_elem_double64* L){

    if(!*isMdiag)
      uni10_linalg::matrixLQ(Mij_ori->__elem, *M, *N, Q->__elem, L->__elem);
    else
      uni10_error_msg(true, "Developping!!!");

  }

  void matrixSVD(const uni10_elem_double64* Mij_ori, uni10_const_bool* isMdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_double64* U, uni10_elem_double64* S, uni10_elem_double64* vT){

    if(!*isMdiag)
      uni10_linalg::matrixSVD(Mij_ori->__elem, *M, *N, U->__elem, S->__elem, vT->__elem);
    else
      uni10_error_msg(true, "Developping!!!");

  }

  void matrixInv(const uni10_elem_double64* A, const uni10_uint64* N, uni10_const_bool* isMdiag){

    if(!*isMdiag)
      uni10_linalg::matrixInv(A->__elem, *N);
    else
      uni10_error_msg(true, "Developping!!!");

  }

  // LAPACK
  //
  // UNI10_COMPLEX128
  void matrixQR(const uni10_elem_complex128* Mij_ori, uni10_const_bool* isMdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_complex128* Q, uni10_elem_complex128* R){

    if(!*isMdiag)
      uni10_linalg::matrixQR(Mij_ori->__elem, *M, *N, Q->__elem, R->__elem);
    else
      uni10_error_msg(true, "Developping!!!");

  }

  void matrixRQ(const uni10_elem_complex128* Mij_ori, uni10_const_bool* isMdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_complex128* Q, uni10_elem_complex128* R){

    if(!*isMdiag)
      uni10_linalg::matrixRQ(Mij_ori->__elem, *M, *N, Q->__elem, R->__elem);
    else
      uni10_error_msg(true, "Developping!!!");

  }

  void matrixQL(const uni10_elem_complex128* Mij_ori, uni10_const_bool* isMdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_complex128* Q, uni10_elem_complex128* L){

    if(!*isMdiag)
      uni10_linalg::matrixQL(Mij_ori->__elem, *M, *N, Q->__elem, L->__elem);
    else
      uni10_error_msg(true, "Developping!!!");

  }

  void matrixLQ(const uni10_elem_complex128* Mij_ori, uni10_const_bool* isMdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_complex128* Q, uni10_elem_complex128* L){

    if(!*isMdiag)
      uni10_linalg::matrixLQ(Mij_ori->__elem, *M, *N, Q->__elem, L->__elem);
    else
      uni10_error_msg(true, "Developping!!!");

  }

  void matrixSVD(const uni10_elem_complex128* Mij_ori, uni10_const_bool* isMdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_complex128* U, uni10_elem_complex128* S, uni10_elem_complex128* vT){

    if(!*isMdiag)
      uni10_linalg::matrixSVD(Mij_ori->__elem, *M, *N, U->__elem, S->__elem, vT->__elem);
    else
      uni10_error_msg(true, "Developping!!!");

  }

  void matrixInv(const uni10_elem_complex128* A, const uni10_uint64* N, uni10_const_bool* isMdiag){

    if(!*isMdiag)
      uni10_linalg::matrixInv(A->__elem, *N);
    else
      uni10_error_msg(true, "Developping!!!");

  }

}
