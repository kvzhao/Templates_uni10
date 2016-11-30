#ifndef __UNI10_LINALG_API_H__
#define __UNI10_LINALG_API_H__

#include <vector>

#include "uni10/uni10_api/Block.h"
#include "uni10/uni10_api/Matrix.h"
//#include "uni10/uni10_api/linalg_typemix.h"

#if defined(CPU) && defined(LAPACK)
#include "uni10/uni10_elem_linalg.h"
#endif

namespace uni10{

  template<typename uni10_type> 
    Matrix<uni10_type> getDiag( const Block<uni10_type>& A ){
      if(A.diag){
        Matrix<uni10_type> D(A);
        return D;
      }
      else{
        Matrix<uni10_type> D(A.Rnum, A.Cnum, true);
        uni10_error_msg(true, "%s", "Developping!!!\n");
        return D;
      }
    }

  template<typename uni10_type> 
    Matrix<uni10_type> dot( const Block<uni10_type>& A, const Block<uni10_type>& B ){

      uni10_error_msg(A.Cnum != B.Rnum, "%s", "The dimensions of the two matrices do not match for matrix multiplication.");
      //Lack of error msgs.
      //
      Matrix<uni10_type> C(A.Rnum, B.Cnum, A.diag && B.diag );

      matrixDot(&A.elem, &A.diag, &B.elem, &B.diag, &A.Rnum, &B.Cnum, &A.Cnum, &C.elem);

      return C;

    }

  template<typename uni10_type>
    std::vector< Matrix<uni10_type> > qr( const Block<uni10_type>& Mij ){

      uni10_error_msg(Mij.Rnum < Mij.Cnum, "%s", "Cannot perform QR decomposition when Rnum < Cnum. Nothing to do." );

      std::vector<Matrix<uni10_type> > outs;
      outs.push_back(Matrix<uni10_type>(Mij.Rnum, Mij.Cnum));
      outs.push_back(Matrix<uni10_type>(Mij.Cnum, Mij.Cnum));

      matrixQR(&Mij.elem, &Mij.diag, &Mij.Rnum, &Mij.Cnum, &outs[0].elem, &outs[1].elem);

      return outs;

    }

  template<typename uni10_type>
    std::vector< Matrix<uni10_type> > rq( const Block<uni10_type>& Mij ){

      uni10_error_msg(Mij.Rnum > Mij.Cnum, "%s", "Cannot perform RQ decomposition when Rnum > Cnum. Nothing to do." );

      std::vector<Matrix<uni10_type> > outs;
      outs.push_back(Matrix<uni10_type>(Mij.Rnum, Mij.Rnum));
      outs.push_back(Matrix<uni10_type>(Mij.Rnum, Mij.Cnum));

      matrixRQ(&Mij.elem, &Mij.diag, &Mij.Rnum, &Mij.Cnum, &outs[0].elem, &outs[1].elem);

      return outs;

    }

  template<typename uni10_type>
    std::vector< Matrix<uni10_type> > lq( const Block<uni10_type>& Mij ){

      uni10_error_msg(Mij.Rnum > Mij.Cnum, "%s", "Cannot perform LQ decomposition when Rnum > Cnum. Nothing to do." );

      std::vector<Matrix<uni10_type> > outs;
      outs.push_back(Matrix<uni10_type>(Mij.Rnum, Mij.Rnum));
      outs.push_back(Matrix<uni10_type>(Mij.Rnum, Mij.Cnum));

      matrixLQ(&Mij.elem, &Mij.diag, &Mij.Rnum, &Mij.Cnum, &outs[0].elem, &outs[1].elem);

      return outs;

    }

  template<typename uni10_type>
    std::vector< Matrix<uni10_type> > ql( const Block<uni10_type>& Mij ){

      uni10_error_msg(Mij.Rnum < Mij.Cnum, "%s", "Cannot perform QL decomposition when Rnum < Cnum. Nothing to do." );

      std::vector<Matrix<uni10_type> > outs;
      outs.push_back(Matrix<uni10_type>(Mij.Rnum, Mij.Cnum));
      outs.push_back(Matrix<uni10_type>(Mij.Cnum, Mij.Cnum));

      matrixQL(&Mij.elem, &Mij.diag, &Mij.Rnum, &Mij.Cnum, &outs[0].elem, &outs[1].elem);

      return outs;

    }

  template<typename uni10_type>
    std::vector< Matrix<uni10_type> > qdr( const Block<uni10_type>& Mij ){

      uni10_error_msg(Mij.Rnum < Mij.Cnum, "%s", "Cannot perform QDR decomposition when Rnum < Cnum. Nothing to do." );

      std::vector<Matrix<uni10_type> > outs;
      outs.push_back(Matrix<uni10_type>(Mij.Rnum, Mij.Cnum));
      outs.push_back(Matrix<uni10_type>(Mij.Cnum, Mij.Cnum, true));
      outs.push_back(Matrix<uni10_type>(Mij.Cnum, Mij.Cnum));

      matrixQDR(&Mij.elem, &Mij.diag, &Mij.Rnum, &Mij.Cnum, &outs[0].elem, &outs[1].elem, &outs[2].elem);

      return outs;

    }

  template<typename uni10_type>
    std::vector< Matrix<uni10_type> > ldq( const Block<uni10_type>& Mij ){

      uni10_error_msg(Mij.Rnum > Mij.Cnum, "%s", "Cannot perform LDQ decomposition when Rnum > Cnum. Nothing to do." );

      std::vector<Matrix<uni10_type> > outs;
      outs.push_back(Matrix<uni10_type>(Mij.Rnum, Mij.Rnum));
      outs.push_back(Matrix<uni10_type>(Mij.Rnum, Mij.Rnum, true));
      outs.push_back(Matrix<uni10_type>(Mij.Rnum, Mij.Cnum));

      matrixLDQ(&Mij.elem, &Mij.diag, &Mij.Rnum, &Mij.Cnum, &outs[0].elem, &outs[1].elem, &outs[2].elem);

      return outs;

    }

  template<typename uni10_type>
    std::vector< Matrix<uni10_type> > qdr_cpivot( const Block<uni10_type>& Mij ){

      uni10_error_msg(Mij.Rnum != Mij.Cnum, "%s", "Cannot perform QDR decomposition when Rnum == Cnum. Nothing to do." );

      std::vector<Matrix<uni10_type> > outs;
      outs.push_back(Matrix<uni10_type>(Mij.Rnum, Mij.Cnum));
      outs.push_back(Matrix<uni10_type>(Mij.Cnum, Mij.Cnum, true));
      outs.push_back(Matrix<uni10_type>(Mij.Cnum, Mij.Cnum));

      matrixQDRCPIVOT(&Mij.elem, &Mij.diag, &Mij.Rnum, &Mij.Cnum, &outs[0].elem, &outs[1].elem, &outs[2].elem);

      return outs;

    }

  template<typename uni10_type>
    std::vector< Matrix<uni10_type> > svd( const Block<uni10_type>& Mij ){

      std::vector<Matrix<uni10_type> > outs;

      uni10_uint64 min = Mij.Rnum < Mij.Cnum ? Mij.Rnum : Mij.Cnum;      //min = min(Rnum,Cnum)
      //GPU_NOT_READY
      outs.push_back(Matrix<uni10_type>(Mij.Rnum, min));
      outs.push_back(Matrix<uni10_type>(min, min, true));
      outs.push_back(Matrix<uni10_type>(min, Mij.Cnum));

      matrixSVD(&Mij.elem, &Mij.diag, &Mij.Rnum, &Mij.Cnum, &outs[0].elem, &outs[1].elem, &outs[2].elem);

      return outs;

    }

  template<typename uni10_type>
    std::vector< Matrix<uni10_type> > eigh( const Block<uni10_type >& Mij ){

      std::vector< Matrix<uni10_type> > outs;
      outs.push_back(Matrix<uni10_type>(Mij.Rnum, Mij.Cnum, true));
      outs.push_back(Matrix<uni10_type>(Mij.Rnum, Mij.Cnum));
      matrixEigh(&Mij.elem, &Mij.diag, &Mij.Cnum, &outs[0].elem, &outs[1].elem);

      return outs;

    }

  template<typename uni10_type>
    std::vector< Matrix<uni10_complex128> > eig( const Block<uni10_type>& Mij ){

      std::vector< Matrix<uni10_complex128> > outs;
      outs.push_back(Matrix<uni10_complex128>(Mij.Rnum, Mij.Cnum, true));
      outs.push_back(Matrix<uni10_complex128>(Mij.Rnum, Mij.Cnum));
      matrixEig(&Mij.elem, &Mij.diag, &Mij.Cnum, &outs[0].elem, &outs[1].elem);

      return outs;

    }


  template<typename uni10_type>
    Matrix<uni10_type> inverse( const Block<uni10_type>& Mij ){

      Matrix<uni10_type> invM(Mij);

      uni10_error_msg(!(Mij.Rnum == Mij.Cnum), "%s", "Cannot perform inversion on a non-square matrix." );

      matrixInv(&invM.elem, &Mij.Rnum, &Mij.diag);

      return invM;

    }

  template<typename uni10_type>
    uni10_type sum( const Block<uni10_type>& Mij ){

      uni10_int32 inc = 1;
      return vectorSum(&Mij.elem, &Mij.elem.__elemNum, &inc);

    }

  template<typename uni10_type>
    uni10_double64 norm( const Block<uni10_type>& Mij ){

      uni10_int32 inc = 1;
      return vectorNorm(&Mij.elem, &Mij.elem.__elemNum, &inc);

    }

  template<typename uni10_type>
    Matrix<uni10_type> transpose( const Block<uni10_type>& Mij ){

      Matrix<uni10_type> MijT(Mij.Cnum, Mij.Rnum, Mij.diag);
      setTranspose(&Mij.elem, &Mij.Rnum, &Mij.Cnum, &MijT.elem);

      return MijT;

    }

  template<typename uni10_type>
    Matrix<uni10_type> dagger( const Block<uni10_type>& Mij ){

      Matrix<uni10_type> MijD(Mij.Cnum, Mij.Rnum, Mij.diag);
      setDagger(&Mij.elem, &Mij.Rnum, &Mij.Cnum, &MijD.elem);

      return MijD;

    }

  template<typename uni10_type>
    Matrix<uni10_type> conj( const Block<uni10_type>& Mij ){

      Matrix<uni10_type> MijConj(Mij.Rnum, Mij.Cnum, Mij.diag);
      setConjugate(&Mij.elem, &Mij.elem.__elemNum, &MijConj.elem);

      return MijConj;

    }

  template<typename uni10_type>
    uni10_type det( const Block<uni10_type>& Mij ){

      uni10_type res = matrixDet(&Mij.elem, &Mij.Rnum, &Mij.diag);

      return res;

    }

  template<typename uni10_type>
    uni10_type trace( const Block<uni10_type>& Mij ){

      uni10_type res = matrixTrace(&Mij.elem, &Mij.diag, &Mij.Rnum, &Mij.Cnum);

      return res;

    }

  template<typename uni10_type>
    Matrix<uni10_type> exph( uni10_double64 a, const Block<uni10_type>& mat){

      std::vector< Matrix<uni10_type> > rets = eigh( mat );

      Matrix<uni10_type> MT = dagger(rets[1]);

      vectorExp( &a, &rets[0].elem, &rets[0].Rnum );

      return MT * (rets[0] * rets[1]);

    }

}

#endif
