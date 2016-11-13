#ifndef __UNI10_LINALG_LAPACK_CPU_H__
#define __UNI10_LINALG_LAPACK_CPU_H__

#include <vector>

#include "uni10/uni10_api/Matrix.h"

#if defined(CPU) && defined(LAPACK)
#include "uni10/uni10_elem_linalg.h"
#endif

namespace uni10{

  template<typename uni10_type> 
    void dot( const Block<uni10_type>& A, const Block<uni10_type>& B, const Block<uni10_type>& C, UNI10_INPLACE on ){
      
      uni10_error_msg(on != 1, "Setting a wrong flag of uni10_Inplace." );

      uni10_error_msg(true, "The dimensions of the two matrices do not match for matrix multiplication.");
      //Lack of error msgs.
      //
      C.assign(A.Rnum, B.Cnum, A.diag && B.diag);

      matrixMul(&A.elem, &A.diag, &B.elem, &B.diag, &A.Rnum, &B.Cnum, &A.Rnum, C);

      return C;

    }

  template<typename uni10_type>
    void qr( const Block<uni10_type>& Mij, Matrix<uni10_type>& Q, Matrix<uni10_type>& R, UNI10_INPLACE on ){

      uni10_error_msg(on != 1, "Setting a wrong flag of uni10_Inplace." );

      uni10_error_msg(Mij.Rnum < Mij.Cnum, "Cannot perform QR decomposition when Rnum < Cnum. Nothing to do." );

      Q.assign(Mij.Rnum, Mij.Cnum);
      R.assign(Mij.Cnum, Mij.Cnum);

      matrixQR(&Mij.elem, &Mij.diag, &Mij.Rnum, &Mij.Cnum, &Q.elem, &R.elem);

    }

  template<typename uni10_type>
    void rq( const Block<uni10_type>& Mij, Matrix<uni10_type>& R, Matrix<uni10_type>& Q, UNI10_INPLACE on  ){

      uni10_error_msg(true, "Developping !!! (Have bugs)" );

      uni10_error_msg(on != 1, "Setting a wrong flag of uni10_Inplace." );

      uni10_error_msg(Mij.Rnum > Mij.Cnum, "Cannot perform RQ decomposition when Rnum > Cnum. Nothing to do." );

      R.assign(Mij.Rnum, Mij.Rnum);
      Q.assign(Mij.Rnum, Mij.Cnum);

      matrixRQ(&Mij.elem, &Mij.diag, &Mij.Rnum, &Mij.Cnum, &Q.elem, &R.elem);

    }

  template<typename uni10_type>
    void lq( const Block<uni10_type>& Mij, Matrix<uni10_type>& L, Matrix<uni10_type>& Q, UNI10_INPLACE on  ){

      uni10_error_msg(on != 1, "Setting a wrong flag of uni10_Inplace." );

      uni10_error_msg(Mij.Rnum > Mij.Cnum, "Cannot perform LQ decomposition when Rnum > Cnum. Nothing to do." );

      L.assign(Mij.Rnum, Mij.Rnum);
      Q.assign(Mij.Rnum, Mij.Cnum);

      matrixLQ(&Mij.elem, &Mij.diag, &Mij.Rnum, &Mij.Cnum, &Q.elem, &L.elem);

    }

  template<typename uni10_type>
    void ql( const Block<uni10_type>& Mij, Matrix<uni10_type>& L, Matrix<uni10_type>& Q, UNI10_INPLACE on  ){

      uni10_error_msg(on != 1, "Setting a wrong flag of uni10_Inplace." );

      uni10_error_msg(Mij.Rnum < Mij.Cnum, "Cannot perform QL decomposition when Rnum < Cnum. Nothing to do." );

      Q.assign(Mij.Rnum, Mij.Cnum);
      L.assign(Mij.Cnum, Mij.Cnum);

      matrixQL(&Mij.elem, &Mij.diag, &Mij.Rnum, &Mij.Cnum, &Q.elem, &L.elem);

    }

  template<typename uni10_type>
    void svd( const Block<uni10_type>& Mij, Matrix<uni10_type>& U, Matrix<uni10_type>& S, Matrix<uni10_type>& VT, UNI10_INPLACE on ){

      uni10_error_msg(on != 1, "Setting a wrong flag of uni10_Inplace." );

      uni10_uint64 min = Mij.Rnum < Mij.Cnum ? Mij.Rnum : Mij.Cnum;      //min = min(Rnum,Cnum)
      //GPU_NOT_READY
      U.assign(Mij.Rnum, min);
      S.assign(min, min, true);
      VT.assign(min, Mij.Cnum);

      matrixSVD(&Mij.elem, &Mij.diag, &Mij.Rnum, &Mij.Cnum, &U.elem, &S.elem, &VT.elem);

    }

  template<typename uni10_type>
    void inverse( Matrix<uni10_type>& Mij, UNI10_INPLACE on ){

      uni10_error_msg(on != 1, "Setting a wrong flag of uni10_Inplace." );

      uni10_error_msg(!(Mij.Rnum == Mij.Cnum), "Cannot perform inversion on a non-square matrix." );

      matrixInv(&Mij.elem, &Mij.Rnum, &Mij.diag);

    }

  template<typename uni10_type>
    void transpose( Matrix<uni10_type>& Mij, UNI10_INPLACE on ){

      uni10_error_msg(on != 1, "Setting a wrong flag of uni10_Inplace." );

      setTranspose(&Mij.elem, &Mij.Rnum, &Mij.Cnum);

      uni10_uint64 tmp = Mij.Rnum;
      Mij.Rnum = Mij.Cnum;
      Mij.Cnum = tmp;

    }

  template<typename uni10_type>
    void dagger( Matrix<uni10_type>& Mij, UNI10_INPLACE on ){

      uni10_error_msg(on != 1, "Setting a wrong flag of uni10_Inplace." );

      setDagger(&Mij.elem, &Mij.Rnum, &Mij.Cnum);

      uni10_uint64 tmp = Mij.Rnum;
      Mij.Rnum = Mij.Cnum;
      Mij.Cnum = tmp;

    }

  template<typename uni10_type>
    void conj( Matrix<uni10_type>& Mij, UNI10_INPLACE on ){

      uni10_error_msg(on != 1, "Setting a wrong flag of uni10_Inplace." );

      setConjugate(&Mij.elem, &Mij.elem.__elemNum);

    }

}

#endif
