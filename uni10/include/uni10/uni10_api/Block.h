/****************************************************************************
*  @file Block.h
*  @license
 *   Universal Tensor Network Library
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
//#include "uni10/uni10_api/Matrix.h"

namespace uni10{

  enum UNI10_INPLACE{
    INPLACE = 1
  };

  template<typename uni10_type>
    class UniTensor;

  template<typename uni10_type>
    class Matrix;

  template<typename uni10_type>
    class Block;

  template<typename uni10_type>
    std::ostream& operator<< (std::ostream& os, const Block<uni10_type>& _b);

  template<typename uni10_type>
    Matrix<uni10_type> operator+(const Block<uni10_type>& Ma, const Block<uni10_type>& Mb); 

  template<typename uni10_type>
    Matrix<uni10_type> operator-(const Block<uni10_type>& Ma, const Block<uni10_type>& Mb); 

  template<typename uni10_type>
    Matrix<uni10_type> operator*(const Block<uni10_type>& Ma, const Block<uni10_type>& Mb); 

  template<typename uni10_type>
    Matrix<uni10_type> operator*(uni10_type a, const Block<uni10_type>& Mb); 

  template<typename uni10_type>
    Matrix<uni10_type> operator*( const Block<uni10_type>& Mb, uni10_type a); 

  template<typename uni10_type>
    uni10_bool operator==(const Block<uni10_type>& m1, const Block<uni10_type>& m2);

  template<typename uni10_type>
    uni10_bool operator!=(const Block<uni10_type>& m1, const Block<uni10_type>& m2);

  template<typename uni10_type>
    class Block{  

      protected:

        UELEM(uni10_elem, _package, _type)<uni10_type> elem;     // pointer to a real matrix

        uni10_uint64 Rnum;

        uni10_uint64 Cnum;

        uni10_bool diag;

      public:

        // Four except enforce functions which are designed for tensor_tools.
        uni10_uint64& row_enforce(){return Rnum;};

        uni10_uint64& col_enforce(){return Cnum;};

        uni10_bool& diag_enforce(){return diag;};

        UELEM(uni10_elem, _package, _type)<uni10_type>& elem_enforce(){return elem;};

        const UELEM(uni10_elem, _package, _type)<uni10_type>& const_elem_enforce()const{return elem;};

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

        uni10_type* getElem()const;     // --> uni10_elem().getElem()

        uni10_type at(uni10_uint64 i, uni10_uint64 j)const;

        friend std::ostream& operator<< <>(std::ostream& os, const Block& _b); // --> uni10_elem().print_elem()

        friend Matrix<uni10_type> operator- <>(const Block& m2, const Block& m1); // Elem-wise subtraction

        friend Matrix<uni10_type> operator+ <>(const Block& m1, const Block& m2); // Elem-wise addition

        friend Matrix<uni10_type> operator* <>(const Block& Ma, const Block& Mb); // Elem-wise multiplication

        friend Matrix<uni10_type> operator* <>(uni10_type a, const Block& Ma);

        friend Matrix<uni10_type> operator* <>(const Block& Ma, uni10_type a);

        friend uni10_bool operator== <>(const Block& m1, const Block& m2);

        friend uni10_bool operator!= <>(const Block& m1, const Block& m2);

        //UNI10_LINALG_RETURN_VALUE
        template<typename _uni10_type> 
          friend Matrix<_uni10_type> getDiag( const Block<_uni10_type>& A );

        template<typename _uni10_type> 
          friend Matrix<_uni10_type> dot( const Block<_uni10_type>& A, const Block<_uni10_type>& B );

        template<typename _uni10_type>
          friend std::vector< Matrix<_uni10_type> > qr( const Block<_uni10_type>& M );

        template<typename _uni10_type>
          friend std::vector< Matrix<_uni10_type> > rq( const Block<_uni10_type>& M );

        template<typename _uni10_type>
          friend std::vector< Matrix<_uni10_type> > lq( const Block<_uni10_type>& M );

        template<typename _uni10_type>
          friend std::vector< Matrix<_uni10_type> > ql( const Block<_uni10_type>& M );

        template<typename _uni10_type>
          friend std::vector< Matrix<_uni10_type> > qdr( const Block<_uni10_type>& M);

        template<typename _uni10_type>
          friend std::vector< Matrix<_uni10_type> > ldq( const Block<_uni10_type>& M );

        template<typename _uni10_type>
          friend std::vector< Matrix<_uni10_type> > qdr_cpivot( const Block<_uni10_type>& M );

        template<typename _uni10_type>
          friend std::vector< Matrix<_uni10_type> > svd( const Block<_uni10_type>& M );

        template<typename _uni10_type>
          friend std::vector< Matrix<_uni10_type> > eigh( const Block<_uni10_type>& _Mij);

        template<typename _uni10_type>
          friend std::vector< Matrix<uni10_complex128> > eig( const Block<_uni10_type>& _Mij);

        template<typename _uni10_type>
          friend _uni10_type sum( const Block<_uni10_type>& Mij );

        template<typename _uni10_type>
          friend uni10_double64 norm( const Block<_uni10_type>& Mij );

        template<typename _uni10_type>
          friend Matrix<_uni10_type> inverse( const Block<_uni10_type>& Mij );

        template<typename _uni10_type>
          friend Matrix<_uni10_type> transpose( const Block<_uni10_type>& Mij );

        template<typename _uni10_type>
          friend Matrix<_uni10_type> dagger( const Block<_uni10_type>& Mij );

        template<typename _uni10_type>
          friend Matrix<_uni10_type> conj( const Block<_uni10_type>& Mij );

        template<typename _uni10_type>
          friend _uni10_type det( const Block<_uni10_type>& _Mij );

        template<typename _uni10_type>
          friend _uni10_type trace( const Block<_uni10_type>& _Mij );

        template<typename uni10_typ>
          friend class Matrix;

        template<typename uni10_typ>
          friend class UniTensor;

    };

  template<typename uni10_type>
    std::ostream& operator<< (std::ostream& os, const Block<uni10_type>& _b){

      fprintf(stdout, "\n%ld x %ld = %ld [ Real ElemNum: %ld ]", _b.Rnum, _b.Cnum, _b.Rnum*_b.Cnum, _b.elemNum());

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
    Matrix<uni10_type> operator+ (const Block<uni10_type>& m1, const Block<uni10_type>& m2){

      uni10_error_msg(m1.Rnum != m2.Rnum || m1.Cnum != m2.Cnum, "%s", "Lack err msg!!!");

      Matrix<uni10_type> m3(m1.Rnum, m1.Cnum, m1.diag && m2.diag);
      matrixAdd(&m1.elem, &m1.diag, &m2.elem, &m2.diag, &m1.Rnum, &m1.Cnum, &m3.elem);

      return m3;
       
    }

  template<typename uni10_type>
    Matrix<uni10_type> operator- (const Block<uni10_type>& m1, const Block<uni10_type>& m2){

      uni10_error_msg(m1.Rnum != m2.Rnum || m1.Cnum != m2.Cnum, "%s", "Lack err msg!!!");

      Matrix<uni10_type> m3(m1.Rnum, m1.Cnum, m1.diag && m2.diag);
      matrixSub(&m1.elem, &m1.diag, &m2.elem, &m2.diag, &m1.Rnum, &m1.Cnum, &m3.elem);

      return m3;
       
    }

  template<typename uni10_type>
    Matrix<uni10_type> operator* (const Block<uni10_type>& m1, const Block<uni10_type>& m2){

      uni10_error_msg(m1.Rnum != m2.Rnum || m1.Cnum != m2.Cnum, "%s", "Lack err msg!!!");

      Matrix<uni10_type> m3(m1.Rnum, m1.Cnum, m1.diag || m2.diag);
      matrixMul(&m1.elem, &m1.diag, &m2.elem, &m2.diag, &m1.Rnum, &m1.Cnum, &m3.elem);

      return m3;
       
    }

  template<typename uni10_type>
    Matrix<uni10_type> operator* (uni10_type a, const Block<uni10_type>& m1){
       Matrix<uni10_type> m2(m1);
       m2 *= a;
       return m2;
    }

  template<typename uni10_type>
    Matrix<uni10_type> operator* (const Block<uni10_type>& m1, uni10_type a){
       return a * m1;
    }

  template<typename uni10_type>
    uni10_bool operator== (const Block<uni10_type>& m1, const Block<uni10_type>& m2){

      if( (m1.Rnum != m2.Rnum) || (m1.Cnum != m2.Cnum) || (m1.diag != m2.diag) )
        return false;

      //Check real part 
      for(int i = 0; i < (int)m1.elem->__elemNum; i++)
        if(UNI10_REAL(m1.elem.__elem[i] )- UNI10_REAL(m2.elem.__elem[i] )> 10E-12)
          return false;

      if(m1.elem.__uni10_typeid == 2)
        for(int i = 0; i < (int)m1.elem->__elemNum; i++)
          if(UNI10_IMAG(m1.elem.__elem[i]) - UNI10_IMAG(m2.elem.__elem[i]) > 10E-12)
            return false;

      return true; 
    }

};

/*
        uni10_double64 trace()const;
*/

#endif
