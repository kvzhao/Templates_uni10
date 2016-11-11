#include "uni10/uni10_lapack_cpu/uni10_elem_linalg_lapack_cpu.h"

namespace uni10{

  // BLAS
  //
  // UNI10_DOUBLE64
  void matrixMul(const uni10_elem_double64* A, const uni10_elem_double64* B, const uni10_uint64* M, const uni10_uint64* N, const uni10_uint64* K, uni10_elem_double64* C){

    uni10_linalg::matrixMul(A->elem, B->elem, *M, *N, *K, C->elem);

  }

  void vectorScal(const uni10_double64* a, uni10_elem_double64* X, const uni10_uint64* N){

    uni10_linalg::vectorScal(*a, X->elem, *N);

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
  void matrixQR(const uni10_elem_double64* Mij_ori, const uni10_uint64* M, const uni10_uint64* N, uni10_elem_double64* Q, uni10_elem_double64* R){
    
    uni10_linalg::matrixQR(Mij_ori->elem, *M, *N, Q->elem, R->elem);

  }

  void matrixRQ(const uni10_elem_double64* Mij_ori, const uni10_uint64* M, const uni10_uint64* N, uni10_elem_double64* Q, uni10_elem_double64* R){

    uni10_error_msg(true, "Developping!!!");
    uni10_linalg::matrixRQ(Mij_ori->elem, *M, *N, Q->elem, R->elem);

  }

  void matrixQL(const uni10_elem_double64* Mij_ori, const uni10_uint64* M, const uni10_uint64* N, uni10_elem_double64* Q, uni10_elem_double64* L){

    uni10_linalg::matrixQL(Mij_ori->elem, *M, *N, Q->elem, L->elem);

  }

  void matrixLQ(const uni10_elem_double64* Mij_ori, const uni10_uint64* M, const uni10_uint64* N, uni10_elem_double64* Q, uni10_elem_double64* L){

    uni10_linalg::matrixLQ(Mij_ori->elem, *M, *N, Q->elem, L->elem);

  }

  void matrixSVD(const uni10_elem_double64* Mij_ori, const uni10_uint64* M, const uni10_uint64* N, uni10_elem_double64* U, uni10_elem_double64* S, uni10_elem_double64* vT){

    uni10_linalg::matrixSVD(Mij_ori->elem, *M, *N, U->elem, S->elem, vT->elem);

  }

  void matrixInv(const uni10_elem_double64* A, const uni10_uint64* N){

    uni10_linalg::matrixInv(A->elem, *N);

  }

  // LAPACK
  //
  // UNI10_COMPLEX128
  void matrixQR(const uni10_elem_complex128* Mij_ori, const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* Q, uni10_elem_complex128* R){

    uni10_linalg::matrixQR(Mij_ori->elem, *M, *N, Q->elem, R->elem);

  }

  void matrixRQ(const uni10_elem_complex128* Mij_ori, const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* Q, uni10_elem_complex128* R){

    uni10_linalg::matrixRQ(Mij_ori->elem, *M, *N, Q->elem, R->elem);

  }

  void matrixQL(const uni10_elem_complex128* Mij_ori, const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* Q, uni10_elem_complex128* L){

    uni10_linalg::matrixQL(Mij_ori->elem, *M, *N, Q->elem, L->elem);

  }

  void matrixLQ(const uni10_elem_complex128* Mij_ori, const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* Q, uni10_elem_complex128* L){

    uni10_linalg::matrixLQ(Mij_ori->elem, *M, *N, Q->elem, L->elem);

  }

  void matrixSVD(const uni10_elem_complex128* Mij_ori, const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* U, uni10_elem_complex128* S, uni10_elem_complex128* vT){

    uni10_linalg::matrixSVD(Mij_ori->elem, *M, *N, U->elem, S->elem, vT->elem);

  }

  void matrixInv(const uni10_elem_complex128* A, const uni10_uint64* N){

    uni10_linalg::matrixInv(A->elem, *N);

  }


}
