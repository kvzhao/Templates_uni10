#ifndef __UNI10_LINALG_CPU_DZ_H__
#define __UNI10_LINALG_CPU_DZ_H__

#include "uni10/uni10_type.h"

namespace uni10{

  namespace uni10_linalg{

    void vectorAdd(uni10_complex128* Y, uni10_double64* X, uni10_uint64 N);

    void vectorScal(uni10_double64 a, uni10_complex128* X, uni10_uint64 N);

    void vectorExp(uni10_double64 a, uni10_complex128* X, uni10_uint64 N);

    void eigDecompose(uni10_double64* Kij, uni10_int32 N, uni10_complex128* Eig, uni10_complex128 *EigVec);

    void eigSyDecompose(uni10_complex128* Kij, uni10_int32 N, uni10_double64* Eig, uni10_complex128* EigVec);

    void matrixSVD(uni10_complex128* Mij_ori, uni10_int32 M, uni10_int32 N, uni10_complex128* U, uni10_double64 *S, uni10_complex128* vT);

  };/* namespace uni10_linalg */

};/* namespace uni10 */

#endif
