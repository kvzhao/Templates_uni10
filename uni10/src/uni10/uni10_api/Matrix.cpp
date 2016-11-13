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
#include <math.h>

#include "uni10/uni10_error.h"
#include "uni10/uni10_api/Matrix.h"

namespace uni10{

  template <typename uni10_type>
    Matrix<uni10_type>::Matrix(): Block<uni10_type>(){};

  template <typename uni10_type>
    Matrix<uni10_type>::Matrix(uni10_uint64 _Rnum, uni10_uint64 _Cnum, uni10_bool _diag): Block<uni10_type>(_Rnum, _Cnum, _diag){
      init();
    };

  template <typename uni10_type>
    Matrix<uni10_type>::Matrix(Matrix const& _m): Block<uni10_type>(_m.Rnum, _m.Cnum, _m.diag){
      init(_m.elem.__elem);
    };

  // Copy constructor
  template <typename uni10_type>
    Matrix<uni10_type>::Matrix(Block<uni10_type> const& _b): Block<uni10_type>(_b.Rnum, _b.Cnum, _b.diag){
      init(_b.elem.__elem);
    };

  template <typename uni10_type>
    Matrix<uni10_type>::Matrix(const std::string& fname){
      const std::string f = fname;
      uni10_error_msg(true, "Developping");
    };

  template <typename uni10_type>
    Matrix<uni10_type>::~Matrix(){};

  template <typename uni10_type>
    void Matrix<uni10_type>::assign(uni10_uint64 _Rnum, uni10_uint64 _Cnum, uni10_bool _isdiag){
      this->diag = _isdiag;
      this->Rnum = _Rnum;
      this->Cnum = _Cnum;
      init();
    };

  template <typename uni10_type>
    void Matrix<uni10_type>::load(const std::string& fname){
      const std::string f = fname;
      uni10_error_msg(true, "Developping");
    };

  template <typename uni10_type>

    void Matrix<uni10_type>::setElem(const uni10_type* src, bool src_ongpu){

      uni10_error_msg( src_ongpu, " The source pointer is on the device. Please install the MAGMA or CUDAONLY gpu version instead.");

      if(this->elem.__uni10_typeid != UNI10_TYPE_ID(uni10_type)){

        uni10_error_msg( true, " Developping !!!");

      }

      this->elem.setElem(src);

    };

  template <typename uni10_type>
    void Matrix<uni10_type>::setElem(const std::vector<uni10_type>& elem, bool src_ongpu){

      uni10_error_msg( src_ongpu, " The source pointer is on the device. Please install the MAGMA or CUDAONLY gpu version instead.");

      if(this->diag == false && this->Rnum*this->Cnum != elem.size()){
        char err[512];
        sprintf(err, "Number of the input elements is: %ld, and it doesn't match to the size of matrix: %ld", elem.size(), this->Rnum*this->Cnum);
        uni10_error_msg(true, err);
      }

      if(this->diag == true && fmin(this->Rnum, this->Cnum) != elem.size()){
        char err[512];
        sprintf(err, "Number of the input elements is: %ld, and it doesn't match to the size of matrix: %.0f", elem.size(), fmin(this->Rnum,this->Cnum));
        uni10_error_msg(true, err);
      }

      setElem(&elem[0], src_ongpu);

    };

  template <typename uni10_type>
    void Matrix<uni10_type>::uni10_elem_free(){
      uni10_error_msg(true, "Developping");
    };

  template <typename uni10_type>
    void Matrix<uni10_type>::set_elem_null(){
      uni10_error_msg(true, "Developping");
    };

  template <typename uni10_type>
    void Matrix<uni10_type>::init(const uni10_type* elem){

      if(this->diag){
        this->elem.init(1, fmin(this->Rnum, this->Cnum), elem);
      }
      else{
        this->elem.init(this->Rnum, this->Cnum, elem);
      }

    };

  template class Matrix<uni10_double64>;
  template class Matrix<uni10_complex128>;

};  /* namespace uni10 */
