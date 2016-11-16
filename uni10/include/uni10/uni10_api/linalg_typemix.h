#ifndef __UNI10_LINALG_TYPEMIX_H__
#define __UNI10_LINALG_TYPEMIX_H__

#include <vector>

#include "uni10/uni10_api/Block.h"
#include "uni10/uni10_api/Matrix.h"

#if defined(CPU) && defined(LAPACK)
#include "uni10/uni10_elem_linalg.h"
#endif


namespace uni10{

  void RIN_COUT(const Matrix<uni10_double64>& RIN, Matrix<uni10_complex128>& COUT){
     COUT.assign(RIN.row(), RIN.col(), RIN.isDiag());
     for(uni10_uint64 i = 0; i < COUT.elemNum(); i++)
       COUT[i] = RIN.getElem()[i];
  }

  void RIN_COUT(const Matrix<uni10_complex128>& RIN, Matrix<uni10_complex128>& COUT){
    COUT = RIN;
  }

  void CIN_ROUT(const Matrix<uni10_double64>& CIN, Matrix<uni10_double64>& ROUT){
    ROUT = CIN;
  }

  void CIN_ROUT(const Matrix<uni10_complex128>& CIN, Matrix<uni10_double64>& ROUT){
     ROUT.assign(CIN.row(), CIN.col(), CIN.isDiag());
     for(uni10_uint64 i = 0; i < ROUT.elemNum(); i++){
       ROUT[i] = CIN.getElem()[i].real();
       uni10_error_msg(CIN.getElem()[i].imag() > 1E-12, "%s", "The image terms in this matrix are not 0. Can't convert to real type");
     }
  }

  Matrix<uni10_complex128> RreturnC(const Matrix<uni10_double64>& RIN){
    Matrix<uni10_complex128> COUT(RIN.row(), RIN.col(), RIN.isDiag());
    for(uni10_uint64 i = 0; i < COUT.elemNum(); i++)
      COUT[i] = RIN.getElem()[i];
    return COUT;
  }

  Matrix<uni10_complex128> RreturnC(const Matrix<uni10_complex128>& RIN){
    return RIN;
  }

  Matrix<uni10_double64> CreturnR(const Matrix<uni10_double64>& CIN){
    return CIN;
  }

  Matrix<uni10_double64> CreturnR(const Matrix<uni10_complex128>& CIN){
    Matrix<uni10_double64> ROUT(CIN.row(), CIN.col(), CIN.isDiag());
     for(uni10_uint64 i = 0; i < ROUT.elemNum(); i++){
       ROUT[i] = CIN.getElem()[i].real();
       uni10_error_msg(CIN.getElem()[i].imag() > 1E-12, "%s", "The image terms in this matrix are not 0. Can't convert to real type");
     }
     return ROUT;
  }

}

#endif
