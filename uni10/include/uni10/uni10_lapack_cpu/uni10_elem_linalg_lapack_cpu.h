#ifndef __UNI10_ELEM_LINALG_LAPACK_CPU_H__
#define __UNI10_ELEM_LINALG_LAPACK_CPU_H__

#include "uni10/uni10_type.h"
#include "uni10/uni10_lapack_cpu/uni10_elem_lapack_cpu.h"
#include "uni10/uni10_lapack_cpu/uni10_linalg_cpu_d.h"
#include "uni10/uni10_lapack_cpu/uni10_linalg_cpu_dz.h"
#include "uni10/uni10_lapack_cpu/uni10_linalg_cpu_z.h"

namespace uni10{

  void matrixQR(const uni10_elem_lapack_cpu<uni10_double64>& Mij_ori, uni10_int32& M, uni10_int32& N, uni10_elem_lapack_cpu<uni10_double64>& Q, uni10_elem_lapack_cpu<uni10_double64>& R);

  void matrixQR(const uni10_elem_lapack_cpu<uni10_complex128>& Mij_ori, uni10_int32& M, uni10_int32& N, uni10_elem_lapack_cpu<uni10_complex128>& Q, uni10_elem_lapack_cpu<uni10_complex128>& R);

/*
    // Blas 
  template<typename uni10_type>
    void matrixMul(uni10_elem_lapack_cpu<uni10_type>* A, uni10_elem_lapack_cpu<uni10_type>* B, uni10_int32 M, uni10_int32 N, uni10_int32 K, uni10_elem_lapack_cpu<uni10_type>* C){

      matrixMul(A->elem, B->elem, M, N, K, C->elem);

    };

  template<typename uni10_type>
    void vectorAdd(uni10_elem_lapack_cpu<uni10_type>* Y, uni10_elem_lapack_cpu<uni10_type>* X, uni10_uint64 N){

      vectorAdd(Y->elem, X->elem, N);

    };   // Y = Y + X

  template<typename uni10_type>
    void vectorScal(uni10_elem_lapack_cpu<uni10_type> a, uni10_elem_lapack_cpu<uni10_type>* X, uni10_uint64 N);	  // X = a * X

  template<typename uni10_type>
    void vectorMul(uni10_elem_lapack_cpu<uni10_type>* Y, uni10_elem_lapack_cpu<uni10_type>* X, uni10_uint64 N);   // Y = Y * X, element-wise multiplication;

  template<typename uni10_type>
    uni10_elem_lapack_cpu<uni10_type> vectorSum(uni10_elem_lapack_cpu<uni10_type>* X, uni10_uint64 N, uni10_int32 inc);

  template<typename uni10_type>
    uni10_elem_lapack_cpu<uni10_type> vectorNorm(uni10_elem_lapack_cpu<uni10_type>* X, uni10_uint64 N, uni10_int32 inc);

  template<typename uni10_type>
    void vectorExp(uni10_elem_lapack_cpu<uni10_type> a, uni10_elem_lapack_cpu<uni10_type>* X, uni10_uint64 N);

  template<typename uni10_type>
    void diagRowMul(uni10_elem_lapack_cpu<uni10_type>* mat, uni10_elem_lapack_cpu<uni10_type>* diag, uni10_uint64 M, uni10_uint64 N);

  template<typename uni10_type>
    void diagColMul(uni10_elem_lapack_cpu<uni10_type>* mat, uni10_elem_lapack_cpu<uni10_type>* diag, uni10_uint64 M, uni10_uint64 N);

  template<typename uni10_type>
    void setTranspose(uni10_elem_lapack_cpu<uni10_type>* A, uni10_uint64 M, uni10_uint64 N, uni10_elem_lapack_cpu<uni10_type>* AT);

  template<typename uni10_type>
    void setTranspose(uni10_elem_lapack_cpu<uni10_type>* A, uni10_uint64 M, uni10_uint64 N);

  template<typename uni10_type>
    void setCTranspose(uni10_elem_lapack_cpu<uni10_type>* A, uni10_uint64 M, uni10_uint64 N, uni10_elem_lapack_cpu<uni10_type> *AT);

  template<typename uni10_type>
    void setCTranspose(uni10_elem_lapack_cpu<uni10_type>* A, uni10_uint64 M, uni10_uint64 N);

  template<typename uni10_type>
    void setIdentity(uni10_elem_lapack_cpu<uni10_type>* elem, uni10_uint64 M, uni10_uint64 N);

  template<typename uni10_type>
    void reshapeElem(uni10_elem_lapack_cpu<uni10_type>* elem, uni10_uint64* transOffset);

    // Lapack
  template<typename uni10_type>
    void eigSyDecompose(uni10_elem_lapack_cpu<uni10_type>* Kij, uni10_int32 N, uni10_elem_lapack_cpu<uni10_type>* Eig, uni10_elem_lapack_cpu<uni10_type>* EigVec);

  template<typename uni10_type>
    void matrixSVD(uni10_elem_lapack_cpu<uni10_type>* Mij_ori, uni10_int32 M, uni10_int32 N, uni10_elem_lapack_cpu<uni10_type>* U, uni10_elem_lapack_cpu<uni10_type>* S, uni10_elem_lapack_cpu<uni10_type>* vT);

  template<typename uni10_type>
    void matrixInv(uni10_elem_lapack_cpu<uni10_type>* A, uni10_int32 N, bool diag);

  template<typename uni10_type>
    void matrixRQ(uni10_elem_lapack_cpu<uni10_type>* Mij_ori, uni10_int32 M, uni10_int32 N, uni10_elem_lapack_cpu<uni10_type>* Q, uni10_elem_lapack_cpu<uni10_type>* R);

  template<typename uni10_type>
    void matrixQL(uni10_elem_lapack_cpu<uni10_type>* Mij_ori, uni10_int32 M, uni10_int32 N, uni10_elem_lapack_cpu<uni10_type>* Q, uni10_elem_lapack_cpu<uni10_type>* L);

  template<typename uni10_type>
    void matrixLQ(uni10_elem_lapack_cpu<uni10_type>* Mij_ori, uni10_int32 M, uni10_int32 N, uni10_elem_lapack_cpu<uni10_type>* Q, uni10_elem_lapack_cpu<uni10_type>* L);
*/

}

#endif

