#ifndef __UNI10_ELEM_LINALG_LAPACK_CPU_H__
#define __UNI10_ELEM_LINALG_LAPACK_CPU_H__

#include "uni10/uni10_type.h"
#include "uni10/uni10_lapack_cpu/uni10_elem_lapack_cpu.h"
#include "uni10/uni10_lapack_cpu/uni10_linalg_cpu_d.h"
#include "uni10/uni10_lapack_cpu/uni10_linalg_cpu_dz.h"
#include "uni10/uni10_lapack_cpu/uni10_linalg_cpu_z.h"

namespace uni10{

  // Blas 
  //
  // UNI10_DOUBLE64
  void vectorAdd(uni10_elem_double64* Y, const uni10_elem_double64* X, const uni10_uint64* N);

  void vectorSub(uni10_elem_double64* Y, const uni10_elem_double64* X, const uni10_uint64* N);

  void vectorMul(uni10_elem_double64* Y, const uni10_elem_double64* X, const uni10_uint64* N);   // Y = Y * X, element-wise multiplication;

  void vectorScal(uni10_double64* a, uni10_elem_double64* X, uni10_uint64* N);   // X = a * X

  void vectorExp(uni10_double64* a, uni10_elem_double64* X, uni10_uint64* N);

  uni10_double64 vectorSum (const uni10_elem_double64* X, const uni10_uint64* N, uni10_int32* inc);

  uni10_double64 vectorNorm(const uni10_elem_double64* X, const uni10_uint64* N, uni10_int32* inc);

  void matrixAdd(const uni10_elem_double64* A, uni10_const_bool* Aisdag, const uni10_elem_double64* B, uni10_const_bool* Bisdiag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_double64* C);

  void matrixSub(const uni10_elem_double64* A, uni10_const_bool* Aisdag, const uni10_elem_double64* B, uni10_const_bool* Bisdiag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_double64* C);

  void matrixMul(const uni10_elem_double64* A, uni10_const_bool* Aisdag, const uni10_elem_double64* B, uni10_const_bool* Bisdiag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_double64* C);

  void matrixDot(const uni10_elem_double64* A, uni10_const_bool* Aisdag, const uni10_elem_double64* B, uni10_const_bool* Bisdiag, 
      const uni10_uint64* M, const uni10_uint64* N, const uni10_uint64* K, uni10_elem_double64* C);

  void setTranspose(const uni10_elem_double64* A, const uni10_uint64* M, const uni10_uint64* N, uni10_elem_double64* AT);

  void setTranspose(uni10_elem_double64* A, uni10_uint64* M, uni10_uint64* N);

  void setDagger(const uni10_elem_double64* A, const uni10_uint64* M, const uni10_uint64* N, uni10_elem_double64* AT);

  void setDagger(uni10_elem_double64* A, uni10_uint64* M, uni10_uint64* N);

  void setConjugate(const uni10_elem_double64* A, const uni10_uint64* N, uni10_elem_double64* A_conj);

  void setConjugate(uni10_elem_double64* A, uni10_uint64* N);

  // Blas 
  //
  // UNI10_COMPLEX128
  void vectorAdd(uni10_elem_complex128* Y, const uni10_elem_complex128* X, const uni10_uint64* N);

  void vectorSub(uni10_elem_complex128* Y, const uni10_elem_complex128* X, const uni10_uint64* N);

  void vectorMul(uni10_elem_complex128* Y, const uni10_elem_complex128* X, const uni10_uint64* N);   // Y = Y * X, element-wise multiplication;

  void vectorScal(uni10_complex128* a, uni10_elem_complex128* X, uni10_uint64* N);   // X = a * X

  void vectorExp(uni10_complex128* a, uni10_elem_complex128* X, uni10_uint64* N);

  uni10_complex128 vectorSum (const uni10_elem_complex128* X, const uni10_uint64* N, uni10_int32* inc);

  uni10_double64   vectorNorm(const uni10_elem_complex128* X, const uni10_uint64* N, uni10_int32* inc);

  void matrixAdd(const uni10_elem_complex128* A, uni10_const_bool* Aisdag, const uni10_elem_double64* B, uni10_const_bool* Bisdiag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_double64* C);

  void matrixSub(const uni10_elem_complex128* A, uni10_const_bool* Aisdag, const uni10_elem_double64* B, uni10_const_bool* Bisdiag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_double64* C);

  void matrixMul(const uni10_elem_complex128* A, uni10_const_bool* Aisdag, const uni10_elem_complex128* B, uni10_const_bool* Bisdag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* C);

  void matrixDot(const uni10_elem_complex128* A, uni10_const_bool* Aisdag, const uni10_elem_complex128* B, uni10_const_bool* Bisdag, 
      const uni10_uint64* M, const uni10_uint64* N, const uni10_uint64* K, uni10_elem_complex128* C);

  void setTranspose(const uni10_elem_complex128* A, const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* AT);

  void setTranspose(uni10_elem_complex128* A, uni10_uint64* M, uni10_uint64* N);

  void setDagger(const uni10_elem_complex128* A, const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* AT);

  void setDagger(uni10_elem_complex128* A, uni10_uint64* M, uni10_uint64* N);

  void setConjugate(const uni10_elem_complex128* A, const uni10_uint64* N, uni10_elem_complex128* A_conj);

  void setConjugate(uni10_elem_complex128* A, uni10_uint64* N);

  // Blas 
  //
  // MIX
  void vectorAdd(uni10_elem_complex128* Y, const uni10_elem_complex128* X, const uni10_uint64* N);

  void vectorSub(uni10_elem_complex128* Y, const uni10_elem_complex128* X, const uni10_uint64* N);

  void vectorMul(uni10_elem_complex128* Y, const uni10_elem_double64* X, const uni10_uint64* N);   // Y = Y * X, element-wise multiplication;

  void vectorScal(uni10_double64* a, uni10_elem_complex128* X, uni10_uint64* N);   // X = a * X

  void vectorExp(uni10_double64* a, uni10_elem_complex128* X, uni10_uint64* N);

  void matrixAdd(const uni10_elem_double64* A, uni10_const_bool* Aisdag, const uni10_elem_complex128* B, uni10_const_bool* Bisdag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* C);

  void matrixSub(const uni10_elem_double64* A, uni10_const_bool* Aisdag, const uni10_elem_complex128* B, uni10_const_bool* Bisdag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* C);

  void matrixMul(const uni10_elem_double64* A, uni10_const_bool* Aisdag, const uni10_elem_complex128* B, uni10_const_bool* Bisdag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* C);

  void matrixDot(const uni10_elem_double64* A, uni10_const_bool* Aisdag, const uni10_elem_complex128* B, uni10_const_bool* Bisdag, 
      const uni10_uint64* M, const uni10_uint64* N, const uni10_uint64* K, uni10_elem_complex128* C);

  void matrixAdd(const uni10_elem_complex128* A, uni10_const_bool* Aisdag, const uni10_elem_double64* B, uni10_const_bool* Bisdag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* C);

  void matrixSub(const uni10_elem_complex128* A, uni10_const_bool* Aisdag, const uni10_elem_double64* B, uni10_const_bool* Bisdag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* C);

  void matrixMul(const uni10_elem_complex128* A, uni10_const_bool* Aisdag, const uni10_elem_double64* B, uni10_const_bool* Bisdag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* C);

  void matrixDot(const uni10_elem_complex128* A, uni10_const_bool* Aisdag, const uni10_elem_double64* B, uni10_const_bool* Bisdag, 
      const uni10_uint64* M, const uni10_uint64* N, const uni10_uint64* K, uni10_elem_complex128* C);

  // LAPACK
  //
  //UNI10_DOUBLE64
  void matrixQR(const uni10_elem_double64* Mij_ori, uni10_const_bool* isdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_double64* Q, uni10_elem_double64* R);

  void matrixRQ(const uni10_elem_double64* Mij_ori, uni10_const_bool* isdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_double64* Q, uni10_elem_double64* R);

  void matrixQL(const uni10_elem_double64* Mij_ori, uni10_const_bool* isdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_double64* Q, uni10_elem_double64* L);

  void matrixLQ(const uni10_elem_double64* Mij_ori, uni10_const_bool* isdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_double64* Q, uni10_elem_double64* L);

  void matrixSVD(const uni10_elem_double64* Mij_ori, uni10_const_bool* isdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_double64* U, uni10_elem_double64* S, uni10_elem_double64* vT);

  void matrixInv(const uni10_elem_double64* A, const uni10_uint64* N, uni10_const_bool* isdiag);

  // LAPACK
  //
  //UNI10_COMPLEX128
  void matrixQR(const uni10_elem_complex128* Mij_ori, uni10_const_bool* isdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_complex128* Q, uni10_elem_complex128* R);

  void matrixRQ(const uni10_elem_complex128* Mij_ori, uni10_const_bool* isdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_complex128* Q, uni10_elem_complex128* R);

  void matrixQL(const uni10_elem_complex128* Mij_ori, uni10_const_bool* isdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_complex128* Q, uni10_elem_complex128* L);

  void matrixLQ(const uni10_elem_complex128* Mij_ori, uni10_const_bool* isdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_complex128* Q, uni10_elem_complex128* L);

  void matrixSVD(const uni10_elem_complex128* Mij_ori, uni10_const_bool* isdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_complex128* U, uni10_elem_complex128* S, uni10_elem_complex128* vT);

  void matrixInv(const uni10_elem_complex128* A, const uni10_uint64* N, uni10_const_bool* isdiag);

}

#endif

