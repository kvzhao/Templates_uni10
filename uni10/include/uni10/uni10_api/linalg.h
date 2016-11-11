#ifndef __UNI10_LINALG_LAPACK_CPU_H__
#define __UNI10_LINALG_LAPACK_CPU_H__

#include <vector>

#include "uni10/uni10_type.h"
#include "uni10/uni10_error.h"
#include "uni10/uni10_api/Matrix.h"

#if defined(CPU) && defined(LAPACK)
#include "uni10/uni10_elem_linalg.h"
#endif

namespace uni10{

  template<typename uni10_type>
    std::vector< Matrix<uni10_type> > qr( const Block<uni10_type>& M ){

      uni10_error_msg(M.row() < M.col(), "Cannot perform QR decomposition when Rnum < Cnum. Nothing to do." );

      std::vector<Matrix<uni10_type> > outs;
      outs.push_back(Matrix<uni10_type>(M.row(), M.row()));
      outs.push_back(Matrix<uni10_type>(M.col(), M.col()));
      if(!M->diag)
        matrixQR(&M.elem, M.row(), M.col(), &outs[0].elem, &outs[1].elem);

    };

  template<typename uni10_type>
    void vectorScal(Block<uni10_type> a, Block<uni10_type>* X, uni10_uint64 N);	  // X = a * X

  template<typename uni10_type>
    void matrixMul(Block<uni10_type>* A, Block<uni10_type>* B, uni10_int32 M, uni10_int32 N, uni10_int32 K, Block<uni10_type>* C);

  template<typename uni10_type>
    void vectorAdd(Block<uni10_type>* Y, Block<uni10_type>* X, uni10_uint64 N);

  template<typename uni10_type>
    void vectorMul(Block<uni10_type>* Y, Block<uni10_type>* X, uni10_uint64 N);   // Y = Y * X, element-wise multiplication;

  template<typename uni10_type>
    Block<uni10_type> vectorSum(Block<uni10_type>* X, uni10_uint64 N, uni10_int32 inc);

  template<typename uni10_type>
    Block<uni10_type> vectorNorm(Block<uni10_type>* X, uni10_uint64 N, uni10_int32 inc);

  template<typename uni10_type>
    void vectorExp(Block<uni10_type> a, Block<uni10_type>* X, uni10_uint64 N);

  template<typename uni10_type>
    void diagRowMul(Block<uni10_type>* mat, Block<uni10_type>* diag, uni10_uint64 M, uni10_uint64 N);

  template<typename uni10_type>
    void diagColMul(Block<uni10_type>* mat, Block<uni10_type>* diag, uni10_uint64 M, uni10_uint64 N);
    /*Generate a set of row vectors which form a othonormal basis
     *For the incoming matrix "elem", the number of row <= the number of column, M <= N
     */
  template<typename uni10_type>
    void setTranspose(Block<uni10_type>* A, uni10_uint64 M, uni10_uint64 N, Block<uni10_type>* AT);

  template<typename uni10_type>
    void setTranspose(Block<uni10_type>* A, uni10_uint64 M, uni10_uint64 N);

  template<typename uni10_type>
    void setCTranspose(Block<uni10_type>* A, uni10_uint64 M, uni10_uint64 N, Block<uni10_type> *AT);

  template<typename uni10_type>
    void setCTranspose(Block<uni10_type>* A, uni10_uint64 M, uni10_uint64 N);

  template<typename uni10_type>
    void setIdentity(Block<uni10_type>* elem, uni10_uint64 M, uni10_uint64 N);

  template<typename uni10_type>
    void reshapeElem(Block<uni10_type>* elem, uni10_uint64* transOffset);

    // Lapack
  template<typename uni10_type>
    void eigSyDecompose(Block<uni10_type>* Kij, uni10_int32 N, Block<uni10_type>* Eig, Block<uni10_type>* EigVec);

  template<typename uni10_type>
    void matrixSVD(Block<uni10_type>* Mij_ori, uni10_int32 M, uni10_int32 N, Block<uni10_type>* U, Block<uni10_type>* S, Block<uni10_type>* vT);

  template<typename uni10_type>
    void matrixInv(Block<uni10_type>* A, uni10_int32 N, bool diag);

  template<typename uni10_type>
    void matrixRQ(Block<uni10_type>* Mij_ori, uni10_int32 M, uni10_int32 N, Block<uni10_type>* Q, Block<uni10_type>* R);

  template<typename uni10_type>
    void matrixQL(Block<uni10_type>* Mij_ori, uni10_int32 M, uni10_int32 N, Block<uni10_type>* Q, Block<uni10_type>* L);

  template<typename uni10_type>
    void matrixLQ(Block<uni10_type>* Mij_ori, uni10_int32 M, uni10_int32 N, Block<uni10_type>* Q, Block<uni10_type>* L);


}

#endif
