/****************************************************************************
 *  @file Matrix.h
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
 *  @brief Header file for Matrix class
 *  @author Yun-Da Hsieh
 *  @date 2014-05-06
 *  @since 0.1.0
 *
 *****************************************************************************/
#ifndef __UNI10_MATRIX_H__
#define __UNI10_MATRIX_H__

#include "uni10/uni10_api/Block.h"

namespace uni10{

  enum UNI10_INPLACE{
    INPLACE = 1
  };

  template <typename uni10_type>
    class Matrix:public Block<uni10_type> {

      public:

        template<typename _uni10_type> 
          friend void dot( Matrix<_uni10_type>& A, const Block<_uni10_type>& B, UNI10_INPLACE on );

        template<typename _uni10_type>
          friend void qr( const Block<_uni10_type>& Mij, Matrix<_uni10_type>& Q, Matrix<_uni10_type>& R, UNI10_INPLACE on );

        template<typename _uni10_type>
          friend void rq( const Block<_uni10_type>& Mij, Matrix<_uni10_type>& R, Matrix<_uni10_type>& Q, UNI10_INPLACE on  );

        template<typename _uni10_type>
          friend void lq( const Block<_uni10_type>& Mij, Matrix<_uni10_type>& L, Matrix<_uni10_type>& Q, UNI10_INPLACE on  );

        template<typename _uni10_type>
          friend void ql( const Block<_uni10_type>& Mij, Matrix<_uni10_type>& L, Matrix<_uni10_type>& Q, UNI10_INPLACE on  );

        template<typename _uni10_type>
          friend void svd( const Block<_uni10_type>& Mij, Matrix<_uni10_type>& U, Matrix<_uni10_type>& S, Matrix<_uni10_type>& VT, UNI10_INPLACE on );

        template<typename _uni10_type>
          friend void inverse( Matrix<_uni10_type>& Mij, UNI10_INPLACE on );

        template<typename _uni10_type>
          friend void transpose( Matrix<_uni10_type>& Mij, UNI10_INPLACE on );

        template<typename _uni10_type>
          friend void dagger( Matrix<_uni10_type>& Mij, UNI10_INPLACE on );

        template<typename _uni10_type>
          friend void conj( Matrix<_uni10_type>& Mij, UNI10_INPLACE on );

        Matrix& operator=(const Matrix& _m){
          this->Rnum = _m.Rnum;
          this->Cnum = _m.Cnum;
          this->diag = _m.diag;
          init(_m.elem.__elem);
          return *this;
        };

        Matrix<uni10_type>& operator+=(const Matrix<uni10_type>& _m){
          vectorAdd(&this->elem, &_m.elem, &_m.elem.__elemNum);
        };

        Matrix<uni10_type>& operator-=(const Matrix<uni10_type>& _m){
          vectorSub(&this->elem, &_m.elem, &_m.elem.__elemNum);
        };

        Matrix<uni10_type>& operator*=(const Matrix<uni10_type>& _m){  // elem-wise multiplication
          vectorMul(&this->elem, &_m.elem, &_m.elem.__elemNum);
        };

        template<typename _uni10_type> 
          friend void resize( Matrix<_uni10_type>& A , uni10_uint64 row, uni10_uint64 col);

        explicit Matrix();

        explicit Matrix(uni10_uint64 _Rnum, uni10_uint64 _Cnum, bool _diag=false);

        explicit Matrix<uni10_type>(Block<uni10_type> const& _b);

        explicit Matrix(const std::string& fname);

        Matrix(Matrix const& _m);

        ~Matrix();

        void assign(uni10_uint64 _Rnum, uni10_uint64 _Cnum, uni10_bool _isdiag = false);

        void load(const std::string& fname);

        void setElem(const uni10_type* elem, bool src_ongpu = false);

        void setElem(const std::vector<uni10_type>& elem, bool src_ongpu = false);

      private:

        void uni10_elem_free();

        void set_elem_null();

        void init(const uni10_type* elem = NULL);

    };

  template<typename uni10_type> 
    void resize( Matrix<uni10_type>& A , uni10_uint64 row, uni10_uint64 col){

      A.elem.resize(row, col, A.Rnum, A.Cnum, A.diag);

    }

};  /* namespace uni10 */

#endif /* MATRIX_H */
