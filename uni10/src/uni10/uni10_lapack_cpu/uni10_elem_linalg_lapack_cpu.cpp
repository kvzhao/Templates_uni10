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

  void matrixMul(const uni10_elem_double64* A, uni10_const_bool* isAdiag, const uni10_elem_double64* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, uni10_uint64* N, uni10_uint64* K, uni10_elem_double64* C){

    if(!*isAdiag && !*isBdiag)
      uni10_linalg::matrixMul(A->__elem, B->__elem, *M, *N, *K, C->__elem);
    else
      uni10_error_msg(true, "Developping!!!");

  }

  void diagRowMul(uni10_elem_double64* mat, uni10_elem_double64* diag_mat, uni10_uint64* M, uni10_uint64* N){
    
    uni10_linalg::diagRowMul(mat->__elem, diag_mat->__elem, *M, *N);
      
  }

  void diagColMul(uni10_elem_double64* mat, uni10_elem_double64* diag_mat, uni10_uint64* M, uni10_uint64* N){

    uni10_linalg::diagColMul(mat->__elem, diag_mat->__elem, *M, *N);

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

  void matrixMul(const uni10_elem_complex128* A, uni10_const_bool* isAdiag, const uni10_elem_complex128* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, uni10_uint64* N, uni10_uint64* K, uni10_elem_complex128* C){

      uni10_linalg::matrixMul(A->__elem, B->__elem, *M, *N, *K, C->__elem);

  }

  void diagRowMul(uni10_elem_complex128* mat, uni10_elem_complex128* diag_mat, uni10_uint64* M, uni10_uint64* N){
    
    uni10_linalg::diagRowMul(mat->__elem, diag_mat->__elem, *M, *N);
      
  }

  void diagColMul(uni10_elem_complex128* mat, uni10_elem_complex128* diag_mat, uni10_uint64* M, uni10_uint64* N){

    uni10_linalg::diagColMul(mat->__elem, diag_mat->__elem, *M, *N);

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

  void matrixMul(const uni10_elem_double64* A, uni10_const_bool* isAdiag, const uni10_elem_complex128* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, uni10_uint64* N, uni10_uint64* K, uni10_elem_complex128* C){

      uni10_linalg::matrixMul(A->__elem, B->__elem, *M, *N, *K, C->__elem);

  }

  void matrixMul(const uni10_elem_complex128* A, uni10_const_bool* isAdiag, const uni10_elem_double64* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, uni10_uint64* N, uni10_uint64* K, uni10_elem_complex128* C){

    if(!*isAdiag && !*isBdiag)
      uni10_linalg::matrixMul(A->__elem, B->__elem, *M, *N, *K, C->__elem);
    else
      uni10_error_msg(true, "Developping!!!");

  }

  //void vectorMul(uni10_elem_double64* Y, uni10_elem_double64* X, uni10_uint64 N){

  //}

  //uni10_elem_double64 vectorSum(uni10_elem_double64* X, uni10_uint64 N, uni10_int32 inc){

  //}

  //uni10_elem_double64 vectorNorm(uni10_elem_double64* X, uni10_uint64 N, uni10_int32 inc){

  //}

  //void vectorExp(uni10_elem_double64 a, uni10_elem_double64* X, uni10_uint64 N){

  //}

  //void diagRowMul(uni10_elem_double64* mat, uni10_elem_double64* diag, uni10_uint64 M, uni10_uint64 N){

  //}

  //void diagColMul(uni10_elem_double64* mat, uni10_elem_double64* diag, uni10_uint64 M, uni10_uint64 N){

  //}

  //void setTranspose(uni10_elem_double64* A, uni10_uint64 M, uni10_uint64 N, uni10_elem_double64* AT){

  //}

  //void setTranspose(uni10_elem_double64* A, uni10_uint64 M, uni10_uint64 N){

  //}

  //void setCTranspose(uni10_elem_double64* A, uni10_uint64 M, uni10_uint64 N, uni10_elem_double64 *AT){

  //}

  //void setCTranspose(uni10_elem_double64* A, uni10_uint64 M, uni10_uint64 N){

  //}

  //void setIdentity(uni10_elem_double64* elem, uni10_uint64 M, uni10_uint64 N){

  //}

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
