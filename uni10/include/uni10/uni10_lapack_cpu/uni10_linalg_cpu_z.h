#ifndef __UNI10_LINALG_CPU_Z_H__
#define __UNI10_LINALG_CPU_Z_H__

#include "uni10/uni10_type.h"

namespace uni10{

  namespace uni10_linalg{

    // Blas
    void matrixMul(uni10_complex128* A, uni10_complex128* B, uni10_int32 M, uni10_int32 N, uni10_int32 K, uni10_complex128* C);

    void vectorAdd(uni10_complex128* Y, uni10_complex128* X, uni10_uint64 N);// Y = Y + X

    uni10_complex128 vectorSum(uni10_complex128* X, uni10_uint64 N, uni10_int32 inc);

    double vectorNorm(uni10_complex128* X, uni10_uint64 N, uni10_int32 inc);

    void setTranspose(uni10_complex128* A, uni10_uint64 M, uni10_uint64 N, uni10_complex128* AT);

    void setTranspose(uni10_complex128* A, uni10_uint64 M, uni10_uint64 N);

    void setCTranspose(uni10_complex128* A, uni10_uint64 M, uni10_uint64 N, uni10_complex128* AT);

    void setCTranspose(uni10_complex128* A, uni10_uint64 M, uni10_uint64 N);

    void vectorScal(const uni10_complex128& a, uni10_complex128* X, uni10_uint64 N);	// X = a * X

    void vectorMul(uni10_complex128* Y, uni10_complex128* X, uni10_uint64 N); // Y = Y * X, element-wise multiplication;

    void diagRowMul(uni10_complex128* mat, uni10_complex128* diag, uni10_uint64 M, uni10_uint64 N);

    void diagColMul(uni10_complex128* mat, uni10_complex128* diag, uni10_uint64 M, uni10_uint64 N);

    void vectorExp(const uni10_complex128& a, uni10_complex128* X, uni10_uint64 N);

    // Lapack
    void matrixSVD(uni10_complex128* Mij_ori, uni10_int32 M, uni10_int32 N, uni10_complex128* U, uni10_complex128* S, uni10_complex128* vT);

    void matrixInv(uni10_complex128* A, uni10_int32 N, bool diag);

    void eigDecompose(uni10_complex128* Kij, uni10_int32 N, uni10_complex128* Eig, uni10_complex128 *EigVec);

    void setConjugate(uni10_complex128 *A, uni10_uint64 N);

    void setIdentity(uni10_complex128* elem, uni10_uint64 M, uni10_uint64 N);

    void matrixQR(uni10_complex128* Mij_ori, uni10_int32 M, uni10_int32 N, uni10_complex128* Q, uni10_complex128* R);

    void matrixRQ(uni10_complex128* Mij_ori, uni10_int32 M, uni10_int32 N, uni10_complex128* Q, uni10_complex128* R);

    void matrixQL(uni10_complex128* Mij_ori, uni10_int32 M, uni10_int32 N, uni10_complex128* Q, uni10_complex128* L);

    void matrixLQ(uni10_complex128* Mij_ori, uni10_int32 M, uni10_int32 N, uni10_complex128* Q, uni10_complex128* L);

  };/* namespace uni10_linalg */

};/* namespace uni10 */

#endif
