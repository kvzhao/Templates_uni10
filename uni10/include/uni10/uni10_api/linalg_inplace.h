#ifndef __UNI10_LINALG_INPLACE_API_H__
#define __UNI10_LINALG_INPLACE_API_H__

#include <vector>

#include "uni10/uni10_api/Matrix.h"

#if defined(CPU) && defined(LAPACK)
#include "uni10/uni10_elem_linalg.h"
#endif

namespace uni10{

  template<typename Mat> 
    void dots(const Mat& _m){

    }

  template<typename Mat, typename... Args> 
    void dots(Mat& _m1, const Mat& _m2, const Args&... args) {
      if(_m1.elemNum() == 0)
        _m1 = _m2;
      else
        dot(_m1, _m2, INPLACE);
      dots(_m1, args...);
    }

  template<typename _uni10_type> 
    void resize( Matrix<_uni10_type>& A , uni10_uint64 row, uni10_uint64 col){ 

      A.elem.resize(row, col, A.Rnum, A.Cnum, A.diag);

    }

  template<typename uni10_type> 
    void dot( Matrix<uni10_type>& A, const Matrix<uni10_type>& B, UNI10_INPLACE on ){
      
      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace." );

      uni10_error_msg(A.Cnum != B.Rnum, "%s", "The dimensions of the two matrices do not match for matrix multiplication.");
      //Lack of error msgs.
      //
      Matrix<uni10_type> C(A.Rnum, B.Cnum, A.diag && B.diag);

      matrixDot(&A.elem, &A.diag, &B.elem, &B.diag, &A.Rnum, &B.Cnum, &A.Cnum, &C.elem);

      A = C;

    }

  template<typename uni10_type> 
    void dot( const Matrix<uni10_type>& A, const Matrix<uni10_type>& B, Matrix<uni10_type>& C, UNI10_INPLACE on ){
      
      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace." );

      uni10_error_msg(A.Cnum != B.Rnum, "%s", "The dimensions of the two matrices do not match for matrix multiplication.");
      //Lack of error msgs.
      //
      C.assign(A.Rnum, B.Cnum, A.diag && B.diag);

      matrixDot(&A.elem, &A.diag, &B.elem, &B.diag, &A.Rnum, &B.Cnum, &A.Cnum, &C.elem);

    }

  template<typename uni10_type>
    void qr( const Matrix<uni10_type>& Mij, Matrix<uni10_type>& Q, Matrix<uni10_type>& R, UNI10_INPLACE on ){

      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace." );

      uni10_error_msg(Mij.Rnum < Mij.Cnum, "%s", "Cannot perform QR decomposition when Rnum < Cnum. Nothing to do." );

      Q.assign(Mij.Rnum, Mij.Cnum);
      R.assign(Mij.Cnum, Mij.Cnum);

      matrixQR(&Mij.elem, &Mij.diag, &Mij.Rnum, &Mij.Cnum, &Q.elem, &R.elem);

    }

  template<typename uni10_type>
    void rq( const Matrix<uni10_type>& Mij, Matrix<uni10_type>& R, Matrix<uni10_type>& Q, UNI10_INPLACE on  ){

      uni10_error_msg(true, "%s", "Developping !!! (Have bugs)" );

      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace." );

      uni10_error_msg(Mij.Rnum > Mij.Cnum, "%s", "Cannot perform RQ decomposition when Rnum > Cnum. Nothing to do." );

      R.assign(Mij.Rnum, Mij.Rnum);
      Q.assign(Mij.Rnum, Mij.Cnum);

      matrixRQ(&Mij.elem, &Mij.diag, &Mij.Rnum, &Mij.Cnum, &R.elem, &Q.elem);

    }

  template<typename uni10_type>
    void lq( const Matrix<uni10_type>& Mij, Matrix<uni10_type>& L, Matrix<uni10_type>& Q, UNI10_INPLACE on  ){

      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace." );

      uni10_error_msg(Mij.Rnum > Mij.Cnum, "%s", "Cannot perform LQ decomposition when Rnum > Cnum. Nothing to do." );

      L.assign(Mij.Rnum, Mij.Rnum);
      Q.assign(Mij.Rnum, Mij.Cnum);

      matrixLQ(&Mij.elem, &Mij.diag, &Mij.Rnum, &Mij.Cnum, &L.elem, &Q.elem);

    }

  template<typename uni10_type>
    void ql( const Matrix<uni10_type>& Mij, Matrix<uni10_type>& Q, Matrix<uni10_type>& L, UNI10_INPLACE on  ){

      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace." );

      uni10_error_msg(Mij.Rnum < Mij.Cnum, "%s", "Cannot perform QL decomposition when Rnum < Cnum. Nothing to do." );

      Q.assign(Mij.Rnum, Mij.Cnum);
      L.assign(Mij.Cnum, Mij.Cnum);

      matrixQL(&Mij.elem, &Mij.diag, &Mij.Rnum, &Mij.Cnum, &Q.elem, &L.elem);

    }

  template<typename uni10_type>
    void qdr( const Matrix<uni10_type>& Mij, Matrix<uni10_type>& Q, Matrix<uni10_type>& D, Matrix<uni10_type>& R, UNI10_INPLACE on ){

      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace." );

      uni10_error_msg(Mij.Rnum < Mij.Cnum, "%s", "Cannot perform QDR decomposition when Rnum < Cnum. Nothing to do." );

      Q.assign(Mij.Rnum, Mij.Cnum);
      D.assign(Mij.Cnum, Mij.Cnum, true);
      R.assign(Mij.Cnum, Mij.Cnum);

      matrixQDR(&Mij.elem, &Mij.diag, &Mij.Rnum, &Mij.Cnum, &Q.elem, &D.elem, &R.elem);

    }

  template<typename uni10_type>
    void ldq( const Matrix<uni10_type>& Mij, Matrix<uni10_type>& L, Matrix<uni10_type>& D, Matrix<uni10_type>& Q, UNI10_INPLACE on  ){

      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace." );

      uni10_error_msg(Mij.Rnum > Mij.Cnum, "%s", "Cannot perform LDQ decomposition when Rnum > Cnum. Nothing to do." );

      L.assign(Mij.Rnum, Mij.Rnum);
      D.assign(Mij.Rnum, Mij.Rnum, true);
      Q.assign(Mij.Rnum, Mij.Cnum);

      matrixLDQ(&Mij.elem, &Mij.diag, &Mij.Rnum, &Mij.Cnum, &L.elem, &D.elem, &Q.elem);

    }

  template<typename uni10_type>
    void qdr_cpivot( const Matrix<uni10_type>& Mij, Matrix<uni10_type>& Q, Matrix<uni10_type>& D, Matrix<uni10_type>& R, UNI10_INPLACE on ){

      uni10_error_msg(true, "%s", "Developping" );
      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace." );
      uni10_error_msg(Mij.Rnum != Mij.Cnum, "%s", "Cannot perform QR decomposition when Rnum != Cnum. Nothing to do." );

      Q.assign(Mij.Rnum, Mij.Cnum);
      D.assign(Mij.Cnum, Mij.Cnum, true);
      R.assign(Mij.Cnum, Mij.Cnum);

      matrixQDRCPIVOT(&Mij.elem, &Mij.diag, &Mij.Rnum, &Mij.Cnum, &Q.elem, &D.elem, &R.elem);

    }



  template<typename uni10_type>
    void svd( const Matrix<uni10_type>& Mij, Matrix<uni10_type>& U, Matrix<uni10_type>& S, Matrix<uni10_type>& VT, UNI10_INPLACE on ){

      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace." );

      uni10_uint64 min = Mij.Rnum < Mij.Cnum ? Mij.Rnum : Mij.Cnum;      //min = min(Rnum,Cnum)
      //GPU_NOT_READY
      U.assign(Mij.Rnum, min);
      S.assign(min, min, true);
      VT.assign(min, Mij.Cnum);

      matrixSVD(&Mij.elem, &Mij.diag, &Mij.Rnum, &Mij.Cnum, &U.elem, &S.elem, &VT.elem);

    }

  template<typename uni10_type>
    void eigh( const Matrix<uni10_type>& Mij, Matrix<uni10_double64>& Eig, Matrix<uni10_type>& EigVec, UNI10_INPLACE on ){

      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace." );
      uni10_error_msg(Mij.Rnum != Mij.Cnum, "%s", "Setting a wrong flag of uni10_Inplace." );

      //GPU_NOT_READY
      Eig.assign(Mij.Rnum, Mij.Cnum, true);
      EigVec.assign(Mij.Rnum, Mij.Cnum);

      matrixEigh(&Mij.elem, &Mij.diag, &Mij.Cnum, &Eig.elem, &EigVec.elem);

    }

  template<typename uni10_type>
    void eig( const Matrix<uni10_type>& Mij, Matrix<uni10_complex128>& Eig, Matrix<uni10_complex128>& EigVec, UNI10_INPLACE on ){

      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace." );
      uni10_error_msg(Mij.Rnum != Mij.Cnum, "%s", "Setting a wrong flag of uni10_Inplace." );

      //GPU_NOT_READY
      Eig.assign(Mij.Rnum, Mij.Cnum, true);
      EigVec.assign(Mij.Rnum, Mij.Cnum);

      matrixEig(&Mij.elem, &Mij.diag, &Mij.Cnum, &Eig.elem, &EigVec.elem);

    }

  template<typename uni10_type>
    void inverse( Matrix<uni10_type>& Mij, UNI10_INPLACE on ){

      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace." );

      uni10_error_msg(!(Mij.Rnum == Mij.Cnum), "%s", "Cannot perform inversion on a non-square matrix." );

      matrixInv(&Mij.elem, &Mij.Rnum, &Mij.diag);

    }

  template<typename uni10_type>
    void transpose( Matrix<uni10_type>& Mij, UNI10_INPLACE on ){

      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace." );

      setTranspose(&Mij.elem, &Mij.Rnum, &Mij.Cnum);

      uni10_uint64 tmp = Mij.Rnum;
      Mij.Rnum = Mij.Cnum;
      Mij.Cnum = tmp;

    }

  template<typename uni10_type>
    void dagger( Matrix<uni10_type>& Mij, UNI10_INPLACE on ){

      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace." );

      setDagger(&Mij.elem, &Mij.Rnum, &Mij.Cnum);

      uni10_uint64 tmp = Mij.Rnum;
      Mij.Rnum = Mij.Cnum;
      Mij.Cnum = tmp;

    }

  template<typename uni10_type>
    void conj( Matrix<uni10_type>& Mij, UNI10_INPLACE on ){

      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace." );

      setConjugate(&Mij.elem, &Mij.elem.__elemNum);

    }

}

#endif
