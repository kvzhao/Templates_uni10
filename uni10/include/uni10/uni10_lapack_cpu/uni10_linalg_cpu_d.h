#ifndef __UNI10_LINALG_CPU_D_H__
#define __UNI10_LINALG_CPU_D_H__

#include "uni10/uni10_type.h"

namespace uni10{

  namespace uni10_linalg{

    // Blas 
    //
    void vectorAdd(uni10_double64* Y, uni10_double64* X, uni10_uint64 N);   // Y = Y + X

    void vectorSub(uni10_double64* Y, uni10_double64* X, uni10_uint64 N);   // Y = Y - X

    void vectorMul(uni10_double64* Y, uni10_double64* X, uni10_uint64 N);   // Y = Y * X, element-wise multiplication;

    void vectorScal(uni10_double64 a, uni10_double64* X, uni10_uint64 N);  // X = a * X

    void vectorExp(uni10_double64 a, uni10_double64* X, uni10_uint64 N);

    uni10_double64 vectorSum(uni10_double64* X, uni10_uint64 N, uni10_int32 inc);

    uni10_double64 vectorNorm(uni10_double64* X, uni10_uint64 N, uni10_int32 inc);

    void matrixDot(uni10_double64* A, uni10_double64* B, uni10_int32 M, uni10_int32 N, uni10_int32 K, uni10_double64* C);

    void diagRowMul(uni10_double64* mat, uni10_double64* diag, uni10_uint64 M, uni10_uint64 N);

    void diagColMul(uni10_double64* mat, uni10_double64* diag, uni10_uint64 M, uni10_uint64 N);

    void setTranspose(uni10_double64* A, uni10_uint64 M, uni10_uint64 N, uni10_double64* AT);

    void setTranspose(uni10_double64* A, uni10_uint64 M, uni10_uint64 N);

    void setDagger(uni10_double64* A, uni10_uint64 M, uni10_uint64 N, uni10_double64 *AT);

    void setDagger(uni10_double64* A, uni10_uint64 M, uni10_uint64 N);

    // Lapack
    //
    /*Generate a set of row vectors which form a othonormal basis
     *For the incoming matrix "elem", the number of row <= the number of column, M <= N
     */
    void matrixSVD(uni10_double64* Mij_ori, uni10_int32 M, uni10_int32 N, uni10_double64* U, uni10_double64* S, uni10_double64* vT);

    void matrixQR(uni10_double64* Mij_ori, uni10_int32 M, uni10_int32 N, uni10_double64* Q, uni10_double64* R);

    void matrixRQ(uni10_double64* Mij_ori, uni10_int32 M, uni10_int32 N, uni10_double64* Q, uni10_double64* R);

    void matrixQL(uni10_double64* Mij_ori, uni10_int32 M, uni10_int32 N, uni10_double64* Q, uni10_double64* L);

    void matrixLQ(uni10_double64* Mij_ori, uni10_int32 M, uni10_int32 N, uni10_double64* Q, uni10_double64* L);

    void matrixInv(uni10_double64* A, uni10_int32 N);

    //=====================================================================================//
    //
    void reshapeElem(uni10_double64* elem, uni10_uint64* transOffset);

    void setIdentity(uni10_double64* elem, uni10_uint64 M, uni10_uint64 N);

    void eigSyDecompose(uni10_double64* Kij, uni10_int32 N, uni10_double64* Eig, uni10_double64* EigVec);

  };/* namespace uni10_linalg */

};/* namespace uni10 */

#endif
