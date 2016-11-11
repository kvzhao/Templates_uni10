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

  template <typename uni10_type>
    class Matrix:public Block<uni10_type> {

      public:

        Matrix& operator=(const Matrix& _m){
          this->Rnum = _m.Rnum;
          this->Cnum = _m.Cnum;
          this->diag = _m.diag;
          init(_m.elem.elem);
          return *this;
        };

        explicit Matrix();

        explicit Matrix(uni10_uint64 _Rnum, uni10_uint64 _Cnum, bool _diag=false);

        explicit Matrix<uni10_type>(Block<uni10_type> const& _b);

        explicit Matrix(const std::string& fname);

        Matrix(Matrix const& _m);

        ~Matrix();

        void assign(uni10_uint64 _Rnum, uni10_uint64 _Cnum);

        void load(const std::string& fname);

        void setElem(const uni10_type* elem, bool src_ongpu = false);

        void setElem(const std::vector<uni10_type>& elem, bool src_ongpu = false);

      private:

        void uni10_elem_free();

        void set_elem_null();

        void init(const uni10_type* elem = NULL);

    };

};  /* namespace uni10 */

#endif /* MATRIX_H */
