/****************************************************************************
*  @file Block.h
*  @license
*    Universal Tensor Network Library
*    Copyright (c) 2013-2014
*    National Taiwan University
*    National Tsing-Hua University

*
*    This file is part of Uni10, the Universal Tensor Network Library.
*
*    Uni10 is free software: you can redistribute it and/or modify
*    it under the terms of the GNU Lesser General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version.
*
*    Uni10 is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU Lesser General Public License for more details.
*
*    You should have received a copy of the GNU Lesser General Public License
*    along with Uni10.  If not, see <http://www.gnu.org/licenses/>.
*  @endlicense
*  @brief Header file for Block class
*  @author Yun-Da Hsieh
*  @author Yun-Hsuan Chou
*  @date 2014-05-06
*  @since 0.1.0
*
*****************************************************************************/
#ifndef __UNI10_BLOCK_H__
#define __UNI10_BLOCK_H__

#include <assert.h>

#include <iostream>
#include <iomanip>
#include <vector>

#include "uni10/uni10_type.h"
#include "uni10/uni10_elem.h"
#include "uni10/uni10_lapack_cpu/uni10_elem_linalg_lapack_cpu.h"

namespace uni10{

  template<typename uni10_type>
    class Matrix;

  template<typename uni10_type>
    class Block;

  template<typename uni10_type>
    std::ostream& operator<< (std::ostream& os, const Block<uni10_type>& _b);
    
  template<typename uni10_type>
    class Block{ public: 

      friend std::ostream& operator<< <>(std::ostream& os, const Block& _b); // --> uni10_elem().print_elem()

      //UNI10_LINALG_RETURN_VALUE
      template<typename _uni10_type>
        friend std::vector< Matrix<_uni10_type> > qr( const Block<_uni10_type>& M );

      template<typename _uni10_type>
        friend std::vector< Matrix<_uni10_type> > rq( const Block<_uni10_type>& M );

      template<typename _uni10_type>
        friend std::vector< Matrix<_uni10_type> > lq( const Block<_uni10_type>& M );

      template<typename _uni10_type>
        friend std::vector< Matrix<_uni10_type> > ql( const Block<_uni10_type>& M );

      template<typename _uni10_type>
        friend std::vector< Matrix<_uni10_type> > svd( const Block<_uni10_type>& M );

      template<typename _uni10_type>
        friend Matrix<_uni10_type> inverse( const Block<_uni10_type>& Mij );

        uni10_double64 operator[](uni10_uint64 idx)const;

        explicit Block();

        explicit Block(uni10_uint64 _Rnum, uni10_uint64 _Cnum, bool _diag = false);

        explicit Block(const Block& _b);

        virtual ~Block();

        uni10_uint64 row()const;

        uni10_uint64 col()const;

        bool isDiag()const;

        void save(const std::string& fname)const;

        bool empty()const;              // --> uni10_elem().empty()

        uni10_uint64 elemNum()const;    // --> uni10_elem().elemNum()

        uni10_int32 typeID()const;      // --> uni10_elem().typeid()

        bool isOngpu()const;            // --> uni10_elem().isOngpu()

        uni10_double64* getElem()const; // --> uni10_elem().getElem()

        uni10_type at(uni10_uint64 i, uni10_uint64 j)const;

        friend class Matrix<uni10_type>;

      protected:

        UELEM(uni10_elem, _package, _type)<uni10_type> elem;     // pointer to a real matrix

        uni10_uint64 Rnum;

        uni10_uint64 Cnum;

        bool diag;

    };

  template<typename uni10_type>
    std::ostream& operator<< (std::ostream& os, const Block<uni10_type>& _b){

      fprintf(stdout, "\n%ld x %ld = %ld", _b.Rnum, _b.Cnum, _b.elemNum());

      if(_b.typeID() == 1)  fprintf(stdout, ", REAL");
      else if(_b.typeID() == 2)   os << ", COMPLEX";

      if(_b.diag)
        fprintf(stdout, ", Diagonal");
      fprintf(stdout, "\n\n");

      _b.elem.print_elem(_b.Rnum, _b.Cnum, _b.diag);

      os << "\n";

      return os;
    }

  template<typename uni10_type>
    std::vector< Matrix<uni10_type> > qr( const Block<uni10_type>& Mij ){

      uni10_error_msg(Mij.Rnum < Mij.Cnum, "Cannot perform QR decomposition when Rnum < Cnum. Nothing to do." );

      std::vector<Matrix<uni10_type> > outs;
      outs.push_back(Matrix<uni10_type>(Mij.Rnum, Mij.Cnum));
      outs.push_back(Matrix<uni10_type>(Mij.Cnum, Mij.Cnum));

      if(!Mij.diag)
        matrixQR(&Mij.elem, &Mij.Rnum, &Mij.Cnum, &outs[0].elem, &outs[1].elem);
      else
        uni10_error_msg(true, "Developping!!!");
      //
      return outs;

    }

  template<typename uni10_type>

    std::vector< Matrix<uni10_type> > rq( const Block<uni10_type>& Mij ){

      uni10_error_msg(Mij.Rnum > Mij.Cnum, "Cannot perform RQ decomposition when Rnum > Cnum. Nothing to do." );

      std::vector<Matrix<uni10_type> > outs;
      outs.push_back(Matrix<uni10_type>(Mij.Rnum, Mij.Rnum));
      outs.push_back(Matrix<uni10_type>(Mij.Rnum, Mij.Cnum));

      if(!Mij.diag)
        matrixRQ(&Mij.elem, &Mij.Rnum, &Mij.Cnum, &outs[1].elem, &outs[0].elem);
      else
        uni10_error_msg(true, "Developping!!!");
      //
      return outs;

    }

  template<typename uni10_type>
    std::vector< Matrix<uni10_type> > lq( const Block<uni10_type>& Mij ){

      uni10_error_msg(Mij.Rnum > Mij.Cnum, "Cannot perform LQ decomposition when Rnum > Cnum. Nothing to do." );

      std::vector<Matrix<uni10_type> > outs;
      outs.push_back(Matrix<uni10_type>(Mij.Rnum, Mij.Rnum));
      outs.push_back(Matrix<uni10_type>(Mij.Rnum, Mij.Cnum));

      if(!Mij.diag)
        matrixLQ(&Mij.elem, &Mij.Rnum, &Mij.Cnum, &outs[1].elem, &outs[0].elem);
      else
        uni10_error_msg(true, "Developping!!!");
      //
      return outs;

    }

  template<typename uni10_type>
    std::vector< Matrix<uni10_type> > ql( const Block<uni10_type>& Mij ){

      uni10_error_msg(Mij.Rnum < Mij.Cnum, "Cannot perform QL decomposition when Rnum < Cnum. Nothing to do." );

      std::vector<Matrix<uni10_type> > outs;
      outs.push_back(Matrix<uni10_type>(Mij.Rnum, Mij.Cnum));
      outs.push_back(Matrix<uni10_type>(Mij.Cnum, Mij.Cnum));

      if(!Mij.diag)
        matrixQL(&Mij.elem, &Mij.Rnum, &Mij.Cnum, &outs[0].elem, &outs[1].elem);
      else
        uni10_error_msg(true, "Developping!!!");
      //
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

      if(!Mij.diag)
        matrixSVD(&Mij.elem, &Mij.Rnum, &Mij.Cnum, &outs[0].elem, &outs[1].elem, &outs[2].elem);
      else
        uni10_error_msg(true, "Developping!!!");
      //
      return outs;

    }

  template<typename uni10_type>
    Matrix<uni10_type> inverse( Block<uni10_type> const& Mij ){

      Matrix<uni10_type> invM(Mij);

      uni10_error_msg(!(Mij.Rnum == Mij.Cnum), "Cannot perform inversion on a non-square matrix." );

      if(!invM.diag)
        matrixInv(&invM.elem, &Mij.Rnum);
      else
        uni10_error_msg(true, "Developping!!!");

      return invM;

    }

//#include "uni10/uni10_api/linalg.h"

};


/*
 *      Move to linalg
 *
        uni10_double64 norm()const;

        Block getDiag()const;

        uni10_double64 trace()const;

        uni10_double64 sum()const;

        std::vector<Block> svd()const;

        std::vector<Block> eig()const;

        std::vector<Block> eigh()const;

        Block inverse()const;

        //friend Matrix operator*(const Block& Ma, const Block& Mb); //R*R C*C R*C C*R

        //friend Matrix operator*(Real a, const Block& Ma);

        //friend Matrix operator*(const Block& Ma, Real a);

        //friend Matrix operator*(const Complex& a, const Block& Ma);

        //friend Matrix operator*(const Block& Ma, const Complex& a);

        //friend Matrix operator+(const Block& Ma, const Block& Mb);

        //friend bool operator==(const Block& m1, const Block& m2);

        //friend bool operator!=(const Block& m1, const Block& m2){return !(m1 == m2);};

        //friend Matrix RDotR(const Block& Ma, const Block& Mb);

        //friend Matrix CDotR(const Block& Ma, const Block& Mb);

        //friend Matrix RDotC(const Block& Ma, const Block& Mb);

        //friend Matrix CDotC(const Block& Ma, const Block& Mb);

*/

#endif
