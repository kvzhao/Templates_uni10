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

#include "uni10/uni10_error.h"
#include "uni10/uni10_api/UniTensor.h"

namespace uni10{

  template <typename uni10_type>
    uni10_int32 UniTensor<uni10_type>::COUNTER = 0;

  template <typename uni10_type>
    uni10_uint64 UniTensor<uni10_type>::ELEMNUM = 0;

  template <typename uni10_type>
    uni10_uint64 UniTensor<uni10_type>::MAXELEMNUM = 0;

  template <typename uni10_type>
    uni10_uint64 UniTensor<uni10_type>::MAXELEMTEN = 0;

  template <typename uni10_type>
    UniTensor<uni10_type>::UniTensor(){};

  template <typename uni10_type>
    UniTensor<uni10_type>::UniTensor(uni10_type val): status(0){ 
      this->init();
      this->U_elem.__elem[0] = val;
    }

  template <typename uni10_type>
    UniTensor<uni10_type>::UniTensor(const std::vector<Bond>& _bonds, const std::string& _name): name(_name), bonds(_bonds), status(0){
      this->init();
    }

  template <typename uni10_type>
    UniTensor<uni10_type>::~UniTensor(){
      ELEMNUM -= U_elemNum;
      COUNTER--;
    }

  template <typename uni10_type>
    const Block<uni10_type>& UniTensor<uni10_type>::const_getBlock()const{
      Qnum q0(0);
      return const_getBlock(q0);
    }

  template <typename uni10_type>
    void UniTensor<uni10_type>::setLabel(const uni10_int32 newLabel, const uni10_uint64 idx){
      uni10_error_msg(labels.size() <= idx, "%s", "The bond index is out  of the range of vector(labels).");
      labels[idx] = newLabel;
    }

  template <typename uni10_type>
    void UniTensor<uni10_type>::setLabel(const std::vector<uni10_int32>& newLabels){
      uni10_error_msg(!(bonds.size() == newLabels.size()), "%s", "The size of input vector(labels) does not match for the number of bonds.");
      labels = newLabels;
    }

  template <typename uni10_type>
    void UniTensor<uni10_type>::setLabel(uni10_int32* newLabels){
      std::vector<int> labels(newLabels, newLabels + bonds.size());
      setLabel(labels);
    }

  template <typename uni10_type>
    const Block<uni10_type>& UniTensor<uni10_type>::const_getBlock(const Qnum& qnum)const{
      typename std::map< Qnum, Block<uni10_type> >::const_iterator it = blocks.find(qnum);
      uni10_error_msg(it == blocks.end(), "%s", "There is no block with the given quantum number ");
      return it->second;
    }

  template <typename uni10_type>
    std::vector<Qnum> UniTensor<uni10_type>::blockQnum()const{
      std::vector<Qnum> keys;
      typename std::map< Qnum, Block<uni10_type> >::const_iterator it = blocks.begin();
      for(; it != blocks.end(); it++)
        keys.push_back(it->first);
      return keys;
    }

  template <typename uni10_type>
    Qnum UniTensor<uni10_type>::blockQnum(uni10_uint64 idx)const{
      uni10_error_msg(!(idx < blocks.size()), "Index exceeds the number of the blocks(%ld).", blocks.size());
      typename std::map<Qnum, Block<uni10_type> >::const_iterator it = blocks.begin();
      for(; it != blocks.end(); it++){
        if(idx == 0)
          return it->first;
        idx--;
      }
      return Qnum(0);
    }

  template <typename uni10_type>
    std::map< Qnum, Matrix<uni10_type> > UniTensor<uni10_type>::getBlocks()const{

      std::map< Qnum, Matrix<uni10_type> > mats;
      typename std::map<Qnum, Block<uni10_type> >::const_iterator it = blocks.begin();
      for(; it != blocks.end(); it++){
        Matrix<uni10_type> mat(it->second.Rnum, it->second.Cnum, it->second.elem.__elem, false);
        mats.insert( std::pair<Qnum, Matrix<uni10_type> >(it->first, mat));
      }
      return mats;

    }

  template <typename uni10_type>
    Matrix<uni10_type> UniTensor<uni10_type>::getBlock(uni10_bool diag)const{
      Qnum q0(0);
      return getBlock(q0, diag);
    }

  template <typename uni10_type>
    Matrix<uni10_type> UniTensor<uni10_type>::getBlock(const Qnum& qnum, uni10_bool diag)const{

      typename std::map< Qnum, Block<uni10_type> >::const_iterator it = blocks.find(qnum);
      uni10_error_msg(it == blocks.end(), "%s", "There is no block with the given quantum number");

      if(diag){
        uni10_error_msg(true, "%s", "Developping!!!");
        //return it->second.getDiag();
      }
      else{
        Matrix<uni10_type> mat(it->second.Rnum, it->second.Cnum, it->second.elem.__elem, false);
        return mat;
      }

    }

  template <typename uni10_type>
    void UniTensor<uni10_type>::init(){

      if(bonds.size()){
        U_elemNum = grouping();
        if(!(blocks.size() > 0)){ //No block in Tensor, Error!
          uni10_error_msg(true, "%s", "There is no symmetry block with the given bonds:\n");
          for(uni10_int32 b = 0; b < (uni10_int32)bonds.size(); b++)
            std::cout<<"    "<<bonds[b];
        }

        labels.assign(bonds.size(), 0);
        for(uni10_int32 b = 0; b < (uni10_int32)bonds.size(); b++)
          labels[b] = b;
        status |= HAVEBOND;

      }
      else{
        Qnum q0(0);
        blocks[q0] = Block<uni10_type>(1, 1);
        RBondNum = 0;
        Rdim = 0;
        Cdim = 0;
        U_elemNum = 1;
        status |= HAVEELEM;
      }

      ELEMNUM += U_elemNum;
      COUNTER++;
      if((uni10_uint64)ELEMNUM > MAXELEMNUM)
        MAXELEMNUM = ELEMNUM;
      if(U_elemNum > MAXELEMTEN)
        MAXELEMTEN = U_elemNum;

      U_elem.init(1, U_elemNum, false);
      initBlocks();

    }

  template <typename uni10_type>
    void UniTensor<uni10_type>::initBlocks(){
      uni10_uint64 offset = 0;
      typename std::map< Qnum, Block<uni10_type> >::iterator it = blocks.begin() ;
      for(; it != blocks.end(); it++ ){
        it->second.elem.__elem = &(this->U_elem.__elem[offset]);
        offset += it->second.Rnum * it->second.Cnum;
      }
    }

  template <typename uni10_type>
    uni10_uint64 UniTensor<uni10_type>::grouping(){

      blocks.clear();
      Qnum q0(0);
      uni10_int32   row_bondNum = 0;
      uni10_int32   col_bondNum = 0;
      Rdim = 1;
      Cdim = 1;
      bool IN_BONDS_BEFORE_OUT_BONDS = true;
      for(uni10_int32 i = 0; i < (uni10_int32)bonds.size(); i++){
        uni10_error_msg(bonds[i].Qdegs.size() > 1 || bonds[i].Qnums.size() > 1, "%s", "Ihis UniTensor has symmetry!!");
        if(bonds[i].type() == BD_IN){
          uni10_error_msg(!(IN_BONDS_BEFORE_OUT_BONDS == true), 
              "%s","Error in the input bond array: BD_OUT bonds must be placed after all BD_IN bonds.");
          Rdim *= bonds[i].Qdegs[0];
          row_bondNum++;
        }
        else{
          Cdim *= bonds[i].Qdegs[0];
          col_bondNum++;
          IN_BONDS_BEFORE_OUT_BONDS = false;
        }
      }
      blocks[q0] = Block<uni10_type>(Rdim, Cdim);

      return Rdim*Cdim;
    }

//std::ostream& operator<< (std::ostream& os, const UniTensor& UniT){
//  try{
//    if(!(UniT.status & UniT.HAVEBOND)){
//      if(UniT.ongpu){
//        if(UniT.typeID() == 1)
//          os<<"\nScalar: " << getElemAt(0, UniT.elem, UniT.ongpu);
//        else if(UniT.typeID() == 2)
//          os<<"\nScalar: " << getElemAt(0, UniT.c_elem, UniT.ongpu);
//        os<<", onGPU";
//      }
//      else{
//        if(UniT.typeID() == 1)
//          os<<"\nScalar: " << UniT.elem[0];
//        else if(UniT.typeID() == 2)
//          os<<"\nScalar: " << UniT.c_elem[0];
//      }
//      os<<"\n\n";
//      return os;
//    }
//    int row = 0;
//    int col = 0;
//    std::vector<Bond>bonds = UniT.bond();
//    for(size_t i = 0; i < bonds.size(); i++)
//      if(bonds[i].type() == BD_IN)
//        row++;
//      else
//        col++;
//    int layer = std::max(row, col);
//    int nmlen = UniT.name.length() + 2;
//    int star = 12 + (14 - nmlen) / 2;
//    for(int s = 0; s < star; s++)
//      os << "*";
//    if(UniT.name.length() > 0)
//      os << " " << UniT.name << " ";
//    for(int s = 0; s < star; s++)
//      os<<"*";
//    os<<std::endl;
//    if(UniT.typeID() == 1)
//      os << "REAL" << std::endl;
//    else if(UniT.typeID() == 2)
//      os << "COMPLEX" << std::endl;
//    if(UniT.ongpu)
//      os<<"\n                 onGPU";
//    os << "\n             ____________\n";
//    os << "            |            |\n";
//    int llab = 0;
//    int rlab = 0;
//    char buf[128];
//    for(int l = 0; l < layer; l++){
//      if(l < row && l < col){
//        llab = UniT.labels[l];
//        rlab = UniT.labels[row + l];
//        sprintf(buf, "    %5d___|%-4d    %4d|___%-5d\n", llab, bonds[l].dim(), bonds[row + l].dim(), rlab);
//        os<<buf;
//      }
//      else if(l < row){
//        llab = UniT.labels[l];
//        sprintf(buf, "    %5d___|%-4d    %4s|\n", llab, bonds[l].dim(), "");
//        os<<buf;
//      }
//      else if(l < col){
//        rlab = UniT.labels[row + l];
//        sprintf(buf, "    %5s   |%4s    %4d|___%-5d\n", "", "", bonds[row + l].dim(), rlab);
//        os << buf;
//      }
//      os << "            |            |   \n";
//    }
//    os << "            |____________|\n";
//
//    os << "\n================BONDS===============\n";
//    for(size_t b = 0; b < bonds.size(); b++){
//      os << bonds[b];
//    }
//    os<<"\n===============BLOCKS===============\n";
//    std::map<Qnum, Block> blocks = UniT.const_getBlocks();
//    bool printElem = true;
//    for (std::map<Qnum, Block>::const_iterator  it = blocks.begin() ; it != blocks.end(); it++ ){
//      os << "--- " << it->first << ": ";// << Rnum << " x " << Cnum << " = " << Rnum * Cnum << " ---\n\n";
//      if((UniT.status & UniT.HAVEELEM) && printElem)
//        os<<it->second;
//      else
//        os<<it->second.row() << " x "<<it->second.col()<<": "<<it->second.elemNum()<<std::endl<<std::endl;
//    }
//    os << "Total elemNum: "<<UniT.m_elemNum<<std::endl;
//    os << "***************** END ****************\n\n";
//  }
//  catch(const std::exception& e){
//    propogate_exception(e, "In function operator<<(std::ostream&, uni10::UniTensor&):");
//  }
//  return os;
//}

  //UniTensor::UniTensor(rflag _tp, const std::vector<Bond>& _bonds, std::vector<int>& _labels, const std::string& _name): name(_name), status(0), bonds(_bonds){
  //  try{
  //    throwTypeError(_tp);
  //    initUniT(_tp);
  //    setLabel(_labels);
  //  }
  //  catch(const std::exception& e){
  //    propogate_exception(e, "In constructor UniTensor::UniTensor(std::vector<Bond>&, std::vector<int>&, std::string& = \"\"):");
  //  }
  //}

  //UniTensor::UniTensor(rflag _tp, const std::vector<Bond>& _bonds, int* _labels, const std::string& _name): name(_name), status(0), bonds(_bonds){
  //  try{
  //    throwTypeError(_tp);
  //    initUniT(_tp);
  //    setLabel(_labels);
  //  }
  //  catch(const std::exception& e){
  //    propogate_exception(e, "In constructor UniTensor::UniTensor(std::vector<Bond>&, int*, std::string& = \"\"):");
  //  }
  //}

  //void UniTensor::setRawElem(const std::vector<Real>& rawElem){
  //  try{
  //    setRawElem(&rawElem[0]);
  //  }
  //  catch(const std::exception& e){
  //    propogate_exception(e, "In function UniTensor::setRawElem(std::vector<double>&):");
  //  }
  //}

  //void UniTensor::setRawElem(rflag tp, const Block& blk){
  //  try{
  //    throwTypeError(tp);
  //    setRawElem(blk.getElem(RTYPE));
  //  }
  //  catch(const std::exception& e){
  //    propogate_exception(e, "In function UniTensor::setRawElem(uni10::rflag, uni10::Block&):");
  //  }
  //}

  //void UniTensor::setRawElem(const Real* rawElem){
  //  try{
  //    if((status & HAVEBOND) == 0){
  //      std::ostringstream err;
  //      err<<"Setting elements to a tensor without bonds is not supported.";
  //      throw std::runtime_error(exception_msg(err.str()));
  //    }

  //    if(typeID() == 2)
  //      this->assign(RTYPE, this->bond());

  //    int bondNum = bonds.size();
  //    std::vector<int> Q_idxs(bondNum, 0);
  //    std::vector<int> Q_Bdims(bondNum, 0);
  //    std::vector<int> sB_idxs(bondNum, 0);
  //    std::vector<int> sB_sBdims(bondNum, 0);
  //    std::vector<int> rAcc(bondNum, 1);
  //    for(int b = 0; b < bondNum; b++)
  //      Q_Bdims[b] = bonds[b].Qnums.size();
  //    for(int b = bondNum - 1; b > 0; b--)
  //      rAcc[b - 1] = rAcc[b] * bonds[b].dim();
  //    int Q_off;
  //    int tmp;
  //    int RQoff, CQoff;
  //    size_t sB_r, sB_c;	//sub-block of a Qidx
  //    size_t sB_rDim, sB_cDim;	//sub-block of a Qidx
  //    size_t B_cDim;
  //    size_t E_off;
  //    int R_off;
  //    Real* work = elem;
  //    if(ongpu){
  //      work = (Real*)malloc(m_elemNum * sizeof(Real));
  //    }
  //    for(std::map<int, size_t>::iterator it = QidxEnc.begin(); it != QidxEnc.end(); it++){
  //      Q_off = it->first;
  //      tmp = Q_off;
  //      for(int b = bondNum - 1; b >= 0; b--){
  //        Q_idxs[b] = tmp % Q_Bdims[b];
  //        tmp /= Q_Bdims[b];
  //      }
  //      R_off = 0;
  //      for(int b = 0; b < bondNum; b++){
  //        R_off += rAcc[b] * bonds[b].offsets[Q_idxs[b]];
  //        sB_sBdims[b] = bonds[b].Qdegs[Q_idxs[b]];
  //      }
  //      RQoff = Q_off / CQdim;
  //      CQoff = Q_off % CQdim;
  //      B_cDim = RQidx2Blk[RQoff]->Cnum;
  //      E_off = (RQidx2Blk[RQoff]->m_elem - elem) + (RQidx2Off[RQoff] * B_cDim) + CQidx2Off[CQoff];
  //      sB_rDim = RQidx2Dim[RQoff];
  //      sB_cDim = CQidx2Dim[CQoff];
  //      sB_idxs.assign(bondNum, 0);
  //      for(sB_r = 0; sB_r < sB_rDim; sB_r++)
  //        for(sB_c = 0; sB_c < sB_cDim; sB_c++){
  //          work[E_off + (sB_r * B_cDim) + sB_c] = rawElem[R_off];
  //          for(int bend = bondNum - 1; bend >= 0; bend--){
  //            sB_idxs[bend]++;
  //            if(sB_idxs[bend] < sB_sBdims[bend]){
  //              R_off += rAcc[bend];
  //              break;
  //            }
  //            else{
  //              R_off -= rAcc[bend] * (sB_idxs[bend] - 1);
  //              sB_idxs[bend] = 0;
  //            }
  //          }
  //        }
  //    }

  //    if(ongpu){
  //      elemCopy(elem, work, m_elemNum * sizeof(Real), ongpu, false);
  //      free(work);
  //    }
  //    status |= HAVEELEM;
  //  }
  //  catch(const std::exception& e){
  //    propogate_exception(e, "In function UniTensor::setRawElem(double*):");
  //  }
  //}

  //void UniTensor::setElem(const Real* _elem, bool _ongpu){
  //  try{
  //    if(typeID() == 2)
  //      this->assign(RTYPE, this->bond());
  //    elemCopy(elem, _elem, m_elemNum * sizeof(Real), ongpu, _ongpu);
  //    status |= HAVEELEM;
  //  }
  //  catch(const std::exception& e){
  //    propogate_exception(e, "In function UniTensor::setElem(double*, bool=false):");
  //  }
  //}

  //void UniTensor::setElem(const std::vector<Real>& _elem, bool _ongpu){
  //  try{
  //    setElem(&_elem[0], _ongpu);
  //  }
  //  catch(const std::exception& e){
  //    propogate_exception(e, "In function UniTensor::setElem(std::vector<double>&, bool=false):");
  //  }
  //}

  //void UniTensor::putBlock(rflag tp, const Block& mat, bool force){

  //  try{

  //    //checkUni10TypeError(tp);
  //    throwTypeError(tp);

  //    Qnum q0(0);

  //    putBlock(RTYPE, q0, mat, force);

  //  }
  //  catch(const std::exception& e){
  //    propogate_exception(e, "In function UniTensor::putBlock(uni10::rflag, uni10::Block&):");
  //  }

  //}

  //void UniTensor::putBlock(rflag tp, const Qnum& qnum, const Block& mat, bool force){

  //  try{

  //    //checkUni10TypeError(tp);
  //    throwTypeError(tp);

  //    std::map<Qnum, Block>::iterator it;

  //    if( !force && mat.typeID() == 2){
  //      std::ostringstream err;
  //      err<<"\n1. Can not put a Complex(CTYPE) Matrix into a Real(RTYPE) UniTensor\n\n2. Or you can turn on the force flag, UniTensor::putBlock(qnum, mat, true) or UniTensor::putBlock(RTYPE, qnum, mat, true). \n";
  //      throw std::runtime_error(exception_msg(err.str()));
  //    }

  //    if(!((it = blocks.find(qnum)) != blocks.end())){
  //      std::ostringstream err;
  //      err<<"There is no block with the given quantum number "<<qnum;
  //      throw std::runtime_error(exception_msg(err.str()));
  //    }

  //    if(force && mat.typeID() == 2){

  //      RtoC(*this);

  //      this->putBlock(CTYPE, qnum, mat);

  //    }
  //    else{

  //      if(!(mat.row() == it->second.Rnum && mat.col() == it->second.Cnum)){
  //        std::ostringstream err;
  //        err<<"The dimension of input matrix does not match for the dimension of the block with quantum number "<<qnum<<std::endl;
  //        err<<"  Hint: Use Matrix::resize(int, int)";
  //        throw std::runtime_error(exception_msg(err.str()));
  //      }

  //      if(mat.m_elem != it->second.m_elem){
  //        if(mat.isDiag()){
  //          elemBzero(it->second.m_elem, it->second.Rnum * it->second.Cnum * sizeof(Real), ongpu);
  //          setDiag(it->second.m_elem, mat.getElem(RTYPE), it->second.Rnum, it->second.Cnum, mat.elemNum(), ongpu, mat.isOngpu());
  //        }
  //        else
  //          elemCopy(it->second.m_elem, mat.getElem(RTYPE), it->second.Rnum * it->second.Cnum * sizeof(Real), ongpu, mat.isOngpu());
  //      }

  //      status |= HAVEELEM;

  //    }

  //  }
  //  catch(const std::exception& e){
  //    propogate_exception(e, "In function UniTensor::putBlock(uni10::rflag, uni10::Qnum&, uni10::Block&):");
  //  }

  //}

  //Matrix UniTensor::getRawElem(rflag tp)const{
  //  try{
  //    throwTypeError(tp);
  //    if(status & HAVEBOND && status & HAVEELEM){
  //      int bondNum = bonds.size();
  //      size_t rowNum = 1;
  //      size_t colNum = 1;
  //      for(std::vector<Bond>::const_iterator it = bonds.begin(); it != bonds.end(); ++it){
  //        if(it->type() == BD_IN)
  //          rowNum *= it->dim();
  //        else
  //          colNum *= it->dim();
  //      }
  //      std::vector<size_t> idxs(bondNum, 0);
  //      int bend;
  //      std::vector<Real> rawElem;
  //      while(1){
  //        rawElem.push_back(at(RTYPE, idxs));
  //        for(bend = bondNum - 1; bend >= 0; bend--){
  //          idxs[bend]++;
  //          if(idxs[bend] < bonds[bend].dim())
  //            break;
  //          else
  //            idxs[bend] = 0;
  //        }
  //        if(bend < 0)
  //          break;
  //      }
  //      return Matrix(rowNum, colNum, &rawElem[0]);
  //    }
  //    else if(status & HAVEELEM)
  //      return Matrix(RTYPE, 1, 1, elem);
  //    else
  //      return Matrix();
  //  }
  //  catch(const std::exception& e){
  //    propogate_exception(e, "In function UniTensor::getRawElem(uni10::rflag ):");
  //    return Matrix();
  //  }
  //}

  //Real* UniTensor::getElem(rflag tp){
  //  try{
  //    throwTypeError(tp);
  //    if(typeID() == 2){
  //      std::ostringstream err;
  //      err<<"This Tensor is COMPLEX. Please use UniTensor::getElem(uni10::cflag ) instead";
  //      throw std::runtime_error(exception_msg(err.str()));
  //    }
  //  }
  //  catch(const std::exception& e){
  //    propogate_exception(e, "In function UniTensor::getElem(uni10::rflag ):");
  //  }
  //  return elem;
  //}

  //std::map<Qnum, Matrix> UniTensor::getBlocks(rflag tp)const{
  //  std::map<Qnum, Matrix> mats;
  //  try{
  //    throwTypeError(tp);
  //    for(std::map<Qnum, Block>::const_iterator it = blocks.begin(); it != blocks.end(); it++){
  //      Matrix mat(it->second.Rnum, it->second.Cnum, it->second.m_elem, false, ongpu);
  //      mats.insert(std::pair<Qnum, Matrix>(it->first, mat));
  //    }
  //  }
  //  catch(const std::exception& e){
  //    propogate_exception(e, "In function UniTensor::getBlocks(uni10::rflag ):");
  //  }
  //  return mats;
  //}

  //Matrix UniTensor::getBlock(rflag tp, bool diag)const{
  //  try{
  //    throwTypeError(tp);
  //    Qnum q0(0);
  //    return getBlock(RTYPE, q0, diag);
  //  }
  //  catch(const std::exception& e){
  //    propogate_exception(e, "In function UniTensor::getBlock(uni10::rflag, bool=false):");
  //    return Matrix();
  //  }
  //}

  //Matrix UniTensor::getBlock(rflag tp, const Qnum& qnum, bool diag)const{
  //  try{
  //    throwTypeError(tp);
  //    std::map<Qnum, Block>::const_iterator it = blocks.find(qnum);
  //    if(it == blocks.end()){
  //      std::ostringstream err;
  //      err<<"There is no block with the given quantum number "<<qnum;
  //      throw std::runtime_error(exception_msg(err.str()));
  //    }
  //    if(diag)
  //      return it->second.getDiag();
  //    else{
  //      Matrix mat(it->second.Rnum, it->second.Cnum, it->second.m_elem, false, ongpu);
  //      return mat;
  //    }
  //  }
  //  catch(const std::exception& e){
  //    propogate_exception(e, "In function UniTensor::getBlock(uni10::rflag, uni10::Qnum&):");
  //    return Matrix(0, 0);
  //  }
  //}

  //void UniTensor::set_zero(rflag tp){
  //  try{
  //    throwTypeError(tp);
  //    elemBzero(elem, m_elemNum * sizeof(Real), ongpu);
  //    status |= HAVEELEM;
  //  }
  //  catch(const std::exception& e){
  //    propogate_exception(e, "In function UniTensor::set_zero(uni10::rflag ):");
  //  }
  //}

  //void UniTensor::set_zero(rflag tp, const Qnum& qnum){
  //  try{
  //    throwTypeError(tp);
  //    std::map<Qnum, Block>::iterator it = blocks.find(qnum);
  //    if(it == blocks.end()){
  //      std::ostringstream err;
  //      err<<"There is no block with the given quantum number "<<qnum;
  //      throw std::runtime_error(exception_msg(err.str()));
  //    }
  //    Block& block = it->second;
  //    elemBzero(block.m_elem, block.Rnum * block.Cnum * sizeof(Real), ongpu);
  //    status |= HAVEELEM;
  //  }
  //  catch(const std::exception& e){
  //    propogate_exception(e, "In function UniTensor::set_zero(uni10::rflag, std::Qnum&):");
  //  }
  //}

  //void UniTensor::identity(rflag tp){
  //  try{
  //    throwTypeError(tp);
  //    std::map<Qnum, Block>::iterator it;
  //    for ( it = blocks.begin() ; it != blocks.end(); it++ )
  //      setIdentity(it->second.m_elem, it->second.Rnum, it->second.Cnum, ongpu);
  //    status |= HAVEELEM;
  //  }
  //  catch(const std::exception& e){
  //    propogate_exception(e, "In function UniTensor::identity(uni10::rflag ):");
  //  }
  //}

  //void UniTensor::identity(rflag tp, const Qnum& qnum){
  //  try{
  //    throwTypeError(tp);
  //    std::map<Qnum, Block>::iterator it = blocks.find(qnum);
  //    if(it == blocks.end()){
  //      std::ostringstream err;
  //      err<<"There is no block with the given quantum number "<<qnum;
  //      throw std::runtime_error(exception_msg(err.str()));
  //    }
  //    Block& block = it->second;
  //    setIdentity(block.m_elem, block.Rnum, block.Cnum, ongpu);
  //    status |= HAVEELEM;
  //  }
  //  catch(const std::exception& e){
  //    propogate_exception(e, "In function UniTensor::identity(uni10::rflag, std::Qnum&):");
  //  }
  //}

  //void UniTensor::randomize(rflag tp){
  //  try{
  //    throwTypeError(tp);
  //    elemRand(elem, m_elemNum, ongpu);
  //    status |= HAVEELEM;
  //  }
  //  catch(const std::exception& e){
  //    propogate_exception(e, "In function UniTensor::randomize(uni10::rflag ):");
  //  }
  //}

  //void UniTensor::orthoRand(rflag tp){
  //  try{
  //    throwTypeError(tp);
  //    std::map<Qnum, Block>::iterator it;
  //    for ( it = blocks.begin() ; it != blocks.end(); it++ )
  //      orthoRandomize(it->second.m_elem, it->second.Rnum, it->second.Cnum, ongpu);
  //    status |= HAVEELEM;
  //  }
  //  catch(const std::exception& e){
  //    propogate_exception(e, "In function UniTensor::orthoRand(uni10::rflag ):");
  //  }
  //}

  //void UniTensor::orthoRand(rflag tp, const Qnum& qnum){
  //  try{
  //    throwTypeError(tp);
  //    std::map<Qnum, Block>::iterator it = blocks.find(qnum);
  //    if(it == blocks.end()){
  //      std::ostringstream err;
  //      err<<"There is no block with the given quantum number "<<qnum;
  //      throw std::runtime_error(exception_msg(err.str()));
  //    }
  //    Block& block = it->second;
  //    orthoRandomize(block.m_elem, block.Rnum, block.Cnum, ongpu);
  //    status |= HAVEELEM;
  //  }
  //  catch(const std::exception& e){
  //    propogate_exception(e, "In function UniTensor::orthoRand(uni10::rflag, std::Qnum&):");
  //  }
  //}

  //UniTensor& UniTensor::transpose(rflag tp){
  //  try{
  //    throwTypeError(tp);
  //    if(!(status & HAVEBOND)){
  //      std::ostringstream err;
  //      err<<"There is no bond in the tensor(scalar) to perform transposition.";
  //      throw std::runtime_error(exception_msg(err.str()));
  //    }
  //    int bondNum = bonds.size();
  //    std::vector<int> rsp_outin(bondNum);
  //    int rbondNum = 0;
  //    for(int b = 0; b < bondNum; b++)
  //      if(bonds[b].type() == BD_IN)
  //        rbondNum++;
  //      else
  //        break;
  //    int cbondNum = bondNum - rbondNum;
  //    for(int b = 0; b < bondNum; b++)
  //      if(b < cbondNum)
  //        rsp_outin[b] = rbondNum + b;
  //      else
  //        rsp_outin[b] = b - cbondNum;
  //    std::vector<int> outLabels(bondNum, 0);
  //    std::vector<Bond> outBonds;
  //    for(size_t b = 0; b < bonds.size(); b++){
  //      outBonds.push_back(bonds[rsp_outin[b]]);
  //      outLabels[b] = labels[rsp_outin[b]];
  //    }
  //    for(int b = 0; b < bondNum; b++){
  //      if(b < cbondNum)
  //        outBonds[b].m_type = BD_IN;
  //      else
  //        outBonds[b].m_type = BD_OUT;
  //    }
  //    UniTensor UniTout(RTYPE, outBonds, name);
  //    UniTout.setLabel(outLabels);
  //    if(status & HAVEELEM){
  //      std::map<Qnum, Block>::iterator it_in;
  //      std::map<Qnum, Block>::iterator it_out;
  //      Real* elem_in;
  //      Real* elem_out;
  //      size_t Rnum, Cnum;
  //      for ( it_in = blocks.begin() ; it_in != blocks.end(); it_in++ ){
  //        it_out = UniTout.blocks.find((it_in->first));
  //        Rnum = it_in->second.Rnum;
  //        Cnum = it_in->second.Cnum;
  //        elem_in = it_in->second.m_elem;
  //        elem_out = it_out->second.m_elem;
  //        setTranspose(elem_in, Rnum, Cnum, elem_out, ongpu, UniTout.ongpu);
  //      }
  //      UniTout.status |= HAVEELEM;
  //    }
  //    *this = UniTout;
  //  }
  //  catch(const std::exception& e){
  //    propogate_exception(e, "In function UniTensor::transpose(uni10::rflag ):");
  //  }
  //  return *this;
  //}

  //UniTensor& UniTensor::permute(rflag tp, int rowBondNum){
  //  try{
  //    throwTypeError(tp);
  //    std::vector<int> ori_labels = labels;
  //    this->permute(RTYPE, ori_labels, rowBondNum);
  //  }
  //  catch(const std::exception& e){
  //    propogate_exception(e, "In function UniTensor::permute(uni10::rflag, int):");
  //  }
  //  return *this;
  //}

  //UniTensor& UniTensor::permute(rflag tp, int* newLabels, int rowBondNum){
  //  try{
  //    throwTypeError(tp);
  //    std::vector<int> _labels(newLabels, newLabels + bonds.size());
  //    this->permute(RTYPE, _labels, rowBondNum);
  //  }
  //  catch(const std::exception& e){
  //    propogate_exception(e, "In function UniTensor::permute(uni10::rflag, int*, int):");
  //  }
  //  return *this;
  //}

  //UniTensor& UniTensor::permute(rflag tp, const std::vector<int>& newLabels, int rowBondNum){
  //  try{
  //    throwTypeError(tp);
  //    if((status & HAVEBOND) == 0){
  //      std::ostringstream err;
  //      err<<"There is no bond in the tensor(scalar) to permute.";
  //      throw std::runtime_error(exception_msg(err.str()));
  //    }
  //    if((labels.size() == newLabels.size()) == 0){
  //      std::ostringstream err;
  //      err<<"The size of the input new labels does not match for the number of bonds.";
  //      throw std::runtime_error(exception_msg(err.str()));
  //    }
  //    int bondNum = bonds.size();
  //    std::vector<int> rsp_outin(bondNum);
  //    int cnt = 0;
  //    for(int i = 0; i < bondNum; i++)
  //      for(int j = 0; j < bondNum; j++)
  //        if(labels[i] == newLabels[j]){
  //          rsp_outin[j] = i;
  //          cnt++;
  //        }
  //    if((cnt == newLabels.size()) == 0){
  //      std::ostringstream err;
  //      err<<"The input new labels do not 1-1 correspond to the labels of the tensor.";
  //      throw std::runtime_error(exception_msg(err.str()));
  //    }
  //    bool inorder = true;
  //    for(int i = 1; i < bondNum; i++)
  //      if(rsp_outin[i] != i){
  //        inorder = false;
  //        break;
  //      }
  //    if(inorder && RBondNum == rowBondNum)	//do nothing
  //      return *this;
  //    else{
  //      std::vector<Bond> outBonds;
  //      bool withoutSymmetry = true;
  //      for(size_t b = 0; b < bonds.size(); b++){
  //        outBonds.push_back(bonds[rsp_outin[b]]);
  //        if(bonds[b].Qnums.size() != 1)
  //          withoutSymmetry = false;
  //      }
  //      for(size_t b = 0; b < bonds.size(); b++){
  //        if(b < rowBondNum)
  //          outBonds[b].change(BD_IN);
  //        else
  //          outBonds[b].change(BD_OUT);
  //      }
  //      UniTensor<uni10_type> UniTout(outBonds, name);
  //      if(status & HAVEELEM){
  //        if(withoutSymmetry){
  //          if(!inorder){
  //            if(ongpu && UniTout.ongpu){
  //              size_t* perInfo = (size_t*)malloc(bondNum * 2 * sizeof(size_t));
  //              std::vector<size_t> newAcc(bondNum);
  //              newAcc[bondNum - 1] = 1;
  //              perInfo[bondNum - 1] = 1;
  //              for(int b = bondNum - 1; b > 0; b--){
  //                newAcc[b - 1] = newAcc[b] * UniTout.bonds[b].Qdegs[0];
  //                perInfo[b - 1] = perInfo[b] * bonds[b].Qdegs[0];
  //              }
  //              for(int b = 0; b < bondNum; b++)
  //                perInfo[bondNum + rsp_outin[b]] = newAcc[b];
  //              Real* des_elem = UniTout.elem;
  //              Real* src_elem = elem;
  //              reshapeElem(src_elem, bondNum, m_elemNum, perInfo, des_elem);
  //              free(perInfo);
  //            }
  //            else{
  //              Real* des_elem = UniTout.elem;
  //              Real* src_elem = elem;
  //              size_t memsize = m_elemNum * sizeof(Real);
  //              if(ongpu){
  //                src_elem = (Real*)elemAllocForce(memsize, false);
  //                elemCopy(src_elem, elem, memsize, false, ongpu);
  //              }
  //              if(UniTout.ongpu)
  //                des_elem = (Real*)elemAllocForce(memsize, false);

  //              std::vector<size_t> transAcc(bondNum);
  //              std::vector<size_t> newAcc(bondNum);
  //              transAcc[bondNum - 1] = 1;
  //              newAcc[bondNum - 1] = 1;
  //              for(int b = bondNum - 1; b > 0; b--)
  //                newAcc[b - 1] = newAcc[b] * UniTout.bonds[b].Qdegs[0];
  //              std::vector<int> bondDims(bondNum);
  //              std::vector<int> idxs(bondNum);
  //              for(int b = 0; b < bondNum; b++){
  //                transAcc[rsp_outin[b]] = newAcc[b];
  //                bondDims[b] = bonds[b].Qdegs[0];
  //                idxs[b] = 0;
  //              }
  //              size_t cnt_ot = 0;
  //              for(size_t i = 0; i < m_elemNum; i++){
  //                des_elem[cnt_ot] = src_elem[i];
  //                for(int bend = bondNum - 1; bend >= 0; bend--){
  //                  idxs[bend]++;
  //                  if(idxs[bend] < bondDims[bend]){
  //                    cnt_ot += transAcc[bend];
  //                    break;
  //                  }
  //                  else{
  //                    cnt_ot -= transAcc[bend] * (idxs[bend] - 1);
  //                    idxs[bend] = 0;
  //                  }
  //                }
  //              }
  //              if(ongpu)
  //                elemFree(src_elem, memsize, false);
  //              if(UniTout.ongpu){
  //                elemCopy(UniTout.elem, des_elem, memsize, UniTout.ongpu, false);
  //                elemFree(des_elem, memsize, false);
  //              }
  //            }
  //          }
  //          else{  //non-symmetry inorder
  //            size_t memsize = m_elemNum * sizeof(Real);
  //            elemCopy(UniTout.elem, elem, memsize, UniTout.ongpu, ongpu);
  //          }
  //        }
  //        UniTout.status |= HAVEELEM;
  //      }
  //      *this = UniTout;
  //      this->setLabel(newLabels);
  //    }
  //  }
  //  catch(const std::exception& e){
  //    propogate_exception(e, "In function UniTensor::permute(uni10::rflag, std::vector<int>&, int):");
  //  }
  //  return *this;
  //}

  //Real UniTensor::at(rflag tp, size_t idx)const{
  //  try{
  //    throwTypeError(tp);
  //    if(!(idx < m_elemNum)){
  //      std::ostringstream err;
  //      err<<"Index exceeds the number of elements("<<m_elemNum<<").";
  //      throw std::runtime_error(exception_msg(err.str()));
  //    }
  //    if(typeID() == 2){
  //      std::ostringstream err;
  //      err<<"This Tensor is COMPLEX. Please use UniTensor::at(cflag, const vector<size_t>&) instead";
  //      throw std::runtime_error(exception_msg(err.str()));
  //    }
  //    return getElemAt(idx, elem, ongpu);
  //  }
  //  catch(const std::exception& e){
  //    propogate_exception(e, "In function UniTensor::at(rflag, size_t):");
  //    return 0;
  //  }
  //}

  //UniTensor& UniTensor::combineBond(rflag tp, const std::vector<int>&cmbLabels){
  //  try{
  //    throwTypeError(tp);
  //    if((status & HAVEBOND) == 0){
  //      std::ostringstream err;
  //      err<<"There is no bond in the tensor to be combined.";
  //      throw std::runtime_error(exception_msg(err.str()));
  //    }
  //    if(!(cmbLabels.size() > 1)){
  //      return *this;
  //    }
  //    std::vector<int> rsp_labels(labels.size(), 0);
  //    std::vector<int> reduced_labels(labels.size() - cmbLabels.size() + 1, 0);

  //    std::vector<int> marked(labels.size(), 0);
  //    std::vector<int> picked(cmbLabels.size(), 0);
  //    for(size_t p = 0; p < cmbLabels.size(); p++){
  //      for(size_t l = 0; l < labels.size(); l++){
  //        if(cmbLabels[p] == labels[l]){
  //          picked[p] = l;
  //          marked[l] = 1;
  //          break;
  //        }
  //      }
  //    }
  //    int mark = 0;
  //    for(size_t m = 0; m < marked.size(); m++)
  //      if(marked[m])
  //        mark++;
  //    if(!(mark == cmbLabels.size())){
  //      std::ostringstream err;
  //      err<<"The input labels do not match for the labels of the tensor.";
  //      throw std::runtime_error(exception_msg(err.str()));
  //    }
  //    int enc = 0;
  //    int enc_r = 0;
  //    std::vector<Bond> newBonds;
  //    int RBnum = 0;
  //    for(size_t l = 0; l < labels.size(); l++){
  //      if(marked[l] && l == picked[0]){
  //        for(size_t ll = 0; ll < cmbLabels.size(); ll++){
  //          rsp_labels[enc] = cmbLabels[ll];
  //          enc++;
  //        }
  //        std::vector<Bond> tmpBonds;
  //        for(size_t p = 0; p < picked.size(); p++)
  //          tmpBonds.push_back(bonds[picked[p]]);
  //        if(bonds[picked[0]].type() == BD_IN)
  //          RBnum += picked.size();
  //        newBonds.push_back(combine(tmpBonds));
  //        reduced_labels[enc_r] = labels[l];
  //        enc_r++;
  //      }
  //      else if(marked[l] == 0){
  //        rsp_labels[enc] = labels[l];
  //        reduced_labels[enc_r] = labels[l];
  //        if(bonds[l].type() == BD_IN)
  //          RBnum++;
  //        newBonds.push_back(bonds[l]);
  //        enc_r++;
  //        enc++;
  //      }
  //    }
  //    this->permute(RTYPE, rsp_labels, RBnum);
  //    UniTensor Tout(RTYPE, newBonds, reduced_labels);

  //    if(status & HAVEELEM)
  //      Tout.setElem(elem, ongpu);

  //    *this = Tout;
  //  }
  //  catch(const std::exception& e){
  //    propogate_exception(e, "In function UniTensor::combineBond(uni10::rflag, std::vector<int>&):");
  //  }
  //  return *this;
  //}

  //void UniTensor::addGate(rflag tp, const std::vector<_Swap>& swaps){
  //  try{
  //    throwTypeError(tp);
  //    if((status & HAVEBOND) == 0){
  //      std::ostringstream err;
  //      err<<"Adding swap gates to a tensor without bonds(scalar).";
  //      throw std::runtime_error(exception_msg(err.str()));
  //    }
  //    if((status & HAVEELEM) == 0){
  //      std::ostringstream err;
  //      err<<"Cannot add swap gates to a tensor before setting its elements.";
  //      throw std::runtime_error(exception_msg(err.str()));
  //    }
  //    int sign = 1;
  //    int bondNum = bonds.size();
  //    std::vector<int> Q_idxs(bondNum, 0);
  //    std::vector<int> Q_Bdims(bondNum, 0);
  //    for(int b = 0; b < bondNum; b++)
  //      Q_Bdims[b] = bonds[b].Qnums.size();
  //    int Q_off;
  //    int tmp;
  //    int RQoff, CQoff;
  //    size_t sB_r, sB_c;	//sub-block of a Qidx
  //    size_t sB_rDim, sB_cDim;	//sub-block of a Qidx
  //    size_t B_cDim;
  //    Real* Eptr;
  //    for(std::map<int, size_t>::iterator it = QidxEnc.begin(); it != QidxEnc.end(); it++){
  //      Q_off = it->first;
  //      tmp = Q_off;
  //      for(int b = bondNum - 1; b >= 0; b--){
  //        Q_idxs[b] = tmp % Q_Bdims[b];
  //        tmp /= Q_Bdims[b];
  //      }
  //      RQoff = Q_off / CQdim;
  //      CQoff = Q_off % CQdim;
  //      B_cDim = RQidx2Blk[RQoff]->Cnum;
  //      Eptr = RQidx2Blk[RQoff]->m_elem + (RQidx2Off[RQoff] * B_cDim) + CQidx2Off[CQoff];
  //      sB_rDim = RQidx2Dim[RQoff];
  //      sB_cDim = CQidx2Dim[CQoff];

  //      int sign01 = 0;
  //      for(size_t i = 0; i < swaps.size(); i++)
  //        sign01 ^= (bonds[swaps[i].b1].Qnums[Q_idxs[swaps[i].b1]].prtF() & bonds[swaps[i].b2].Qnums[Q_idxs[swaps[i].b2]].prtF());
  //      sign = sign01 ? -1 : 1;

  //      for(sB_r = 0; sB_r < sB_rDim; sB_r++)
  //        for(sB_c = 0; sB_c < sB_cDim; sB_c++)
  //          Eptr[(sB_r * B_cDim) + sB_c] *= sign;
  //    }
  //  }
  //  catch(const std::exception& e){
  //    propogate_exception(e, "In function UniTensor::addGate(uni10::rflag, std::vector<_Swap>&):");
  //  }
  //}

  //Real UniTensor::trace(rflag tp)const{
  //  try{
  //    throwTypeError(tp);
  //    if(!(status & HAVEELEM)){
  //      std::ostringstream err;
  //      err<<"Cannot trace a tensor before setting its elements.";
  //      throw std::runtime_error(exception_msg(err.str()));
  //    }
  //    if(status & HAVEBOND){
  //      Real trVal = 0;
  //      for(std::map<Qnum, Block>::const_iterator it = blocks.begin() ; it != blocks.end(); it++ ){
  //        if(!(it->second.Rnum == it->second.Cnum)){
  //          std::ostringstream err;
  //          err<<"Cannot trace a non-square block.";
  //          throw std::runtime_error(exception_msg(err.str()));
  //        }
  //        trVal += it->second.trace(RTYPE);
  //      }
  //      return trVal;
  //    }
  //    else{
  //      return getElemAt(0, elem, ongpu);
  //    }
  //  }
  //  catch(const std::exception& e){
  //    propogate_exception(e, "In function UniTensor::trace(uni10::rflag ):");
  //    return 0;
  //  }
  //}

  //UniTensor& UniTensor::partialTrace(rflag tp, int la, int lb){
  //  try{
  //    throwTypeError(tp);
  //    if(!(status & HAVEELEM)){
  //      std::ostringstream err;
  //      err<<"Cannot trace bonds of a tensor before setting its elements.";
  //      throw std::runtime_error(exception_msg(err.str()));
  //    }
  //    if(!(bonds.size() > 2)){
  //      std::ostringstream err;
  //      err<<"The number of bonds must larger than 2 for performing partialTrace.";
  //      throw std::runtime_error(exception_msg(err.str()));
  //    }
  //    int bondNum = bonds.size();
  //    std::vector<Bond> newBonds;
  //    std::vector<int>newLabels(bondNum - 2, 0);
  //    std::vector<int>rsp_labels(bondNum);
  //    int ia, ib;
  //    int enc = 0;
  //    for(size_t l = 0; l < labels.size(); l++){
  //      if(labels[l] == la)
  //        ia = l;
  //      else if(labels[l] == lb)
  //        ib = l;
  //      else{
  //        newBonds.push_back(bonds[l]);
  //        newLabels[enc] = labels[l];
  //        rsp_labels[enc] = labels[l];
  //        enc++;
  //      }
  //    }
  //    if(!(enc == newLabels.size())){
  //      std::ostringstream err;
  //      err<<"Cannot find the two bonds with the given two labels.";
  //      throw std::runtime_error(exception_msg(err.str()));
  //    }

  //    UniTensor Tt(RTYPE, newBonds, newLabels);
  //    rsp_labels[bondNum - 2] = labels[ia];
  //    rsp_labels[bondNum - 1] = labels[ib];
  //    ia = bondNum - 2;
  //    ib = bondNum - 1;
  //    this->permute(RTYPE, rsp_labels, Tt.RBondNum);
  //    std::vector<int> Q_acc(bondNum, 1);
  //    for(int b = bondNum - 1; b > 0; b--)
  //      Q_acc[b - 1] = Q_acc[b] * bonds[b].Qnums.size();
  //    int tQdim = bonds[ia].Qnums.size();
  //    /*Sanity Check*/
  //    if(  !(tQdim == bonds[ib].Qnums.size())  ){
  //      std::ostringstream err;
  //      err<<"The bonds of the given two labels does not match for trace.";
  //      throw std::runtime_error(exception_msg(err.str()));
  //    }
  //    Qnum q0(0, PRT_EVEN);
  //    for(int q = 0; q < tQdim; q++){
  //      if(!((bonds[ia].Qnums[q] * bonds[ib].Qnums[q] == q0) && (bonds[ia].Qdegs[q] == bonds[ib].Qdegs[q]))){
  //        std::ostringstream err;
  //        err<<"The bonds of the given two labels does not match for trace.";
  //        throw std::runtime_error(exception_msg(err.str()));
  //      }
  //    }
  //    /*END*/
  //    int tBnum = Tt.bonds.size();
  //    std::vector<int> Qt_Bdims(tBnum, 0);
  //    for(int b = 0; b < tBnum; b++)
  //      Qt_Bdims[b] = Tt.bonds[b].Qnums.size();

  //    int Qt_off;
  //    int Q_off;
  //    int Qt_RQoff, Qt_CQoff;
  //    int Q_RQoff, Q_CQoff;
  //    size_t sBt_rDim, sBt_cDim;	//sub-block of a Qidx of Tt
  //    size_t sB_rDim, sB_cDim;	//sub-block of a Qidx
  //    size_t Bt_cDim;
  //    Real* Et_ptr;
  //    std::vector<Real*> E_offs(tQdim);
  //    std::vector<size_t> B_cDims(tQdim);
  //    int tQdim2 = tQdim * tQdim;
  //    int Qenc = Q_acc[ia] + Q_acc[ib];
  //    for(std::map<int, size_t>::iterator it = Tt.QidxEnc.begin(); it != Tt.QidxEnc.end(); it++){
  //      Qt_off = it->first;
  //      Qt_RQoff = Qt_off / Tt.CQdim;
  //      Qt_CQoff = Qt_off % Tt.CQdim;
  //      Bt_cDim = Tt.RQidx2Blk[Qt_RQoff]->Cnum;
  //      Et_ptr = Tt.RQidx2Blk[Qt_RQoff]->m_elem + (Tt.RQidx2Off[Qt_RQoff] * Bt_cDim) + Tt.CQidx2Off[Qt_CQoff];
  //      sBt_rDim = Tt.RQidx2Dim[Qt_RQoff];
  //      sBt_cDim = Tt.CQidx2Dim[Qt_CQoff];

  //      for(int q = 0; q < tQdim; q++){
  //        Q_off = Qt_off * tQdim2 + q * Qenc;
  //        Q_RQoff = Q_off / CQdim;
  //        Q_CQoff = Q_off % CQdim;
  //        B_cDims[q] = RQidx2Blk[Q_RQoff]->Cnum;
  //        E_offs[q] = RQidx2Blk[Q_RQoff]->m_elem + (RQidx2Off[Q_RQoff] * B_cDims[q]) + CQidx2Off[Q_CQoff];
  //      }
  //      int tQdeg, sB_c_off;
  //      Real trVal;
  //      for(size_t sB_r = 0; sB_r < sBt_rDim; sB_r++)
  //        for(size_t sB_c = 0; sB_c < sBt_cDim; sB_c++){
  //          trVal = 0;
  //          for(int q = 0; q < tQdim; q++){
  //            tQdeg = bonds[ia].Qdegs[q];
  //            sB_c_off = sB_c * (tQdeg * tQdeg);
  //            for(int t = 0; t < tQdeg; t++){
  //              trVal += E_offs[q][(sB_r * B_cDims[q]) + sB_c_off + t * (tQdeg + 1)];
  //            }
  //          }
  //          Et_ptr[sB_r * Bt_cDim + sB_c] = trVal;
  //        }
  //      Tt.status |= HAVEELEM;
  //    }
  //    *this = Tt;
  //  }
  //  catch(const std::exception& e){
  //    propogate_exception(e, "In function UniTensor::partialTrace(uni10::rflag, int, int):");
  //  }
  //  return *this;
  //}

  //UniTensor& UniTensor::assign(rflag tp, const std::vector<Bond>& _bond){
  //  try{
  //    throwTypeError(tp);
  //    UniTensor T(RTYPE, _bond);
  //    *this = T;
  //  }
  //  catch(const std::exception& e){
  //    propogate_exception(e, "In function UniTensor::assign(uni10::rflag, std::vector<Bond>&):");
  //  }
  //  return *this;
  //}

  //Real UniTensor::at(rflag tp, const std::vector<int>& idxs)const{
  //  try{
  //    throwTypeError(tp);
  //    std::vector<size_t> _idxs(idxs.size());
  //    for(size_t i = 0; i < idxs.size(); i++)
  //      _idxs[i] = idxs[i];
  //    return at(RTYPE, _idxs);
  //  }
  //  catch(const std::exception& e){
  //    propogate_exception(e, "In function UniTensor::at(uni10::rflag, std::vector<int>&):");
  //    return 0;
  //  }
  //}

  //Real UniTensor::at(rflag tp, const std::vector<size_t>& idxs)const{
  //  try{
  //    throwTypeError(tp);
  //    if((status & HAVEBOND) == 0){
  //      std::ostringstream err;
  //      err<<"The tensor is a scalar. Use UniTensor::operator[] instead.";
  //      throw std::runtime_error(exception_msg(err.str()));
  //    }
  //    if(!(idxs.size() == bonds.size())){
  //      std::ostringstream err;
  //      err<<"The size of input indices array does not match with the number of the bonds.";
  //      throw std::runtime_error(exception_msg(err.str()));
  //    }

  //    if(typeID() == 2){
  //      std::ostringstream err;
  //      err<<"This Tensor is COMPLEX. Please use UniTensor::at(uni10::cflag, const vector<size_t>&) instead";
  //      throw std::runtime_error(exception_msg(err.str()));
  //    }

  //    int bondNum = bonds.size();
  //    std::vector<int> Qidxs(bondNum, 0);
  //    for(int b = 0; b < bondNum; b++){
  //      if(!(idxs[b] < bonds[b].dim())){
  //        std::ostringstream err;
  //        err<<"The input indices are out of range.";
  //        throw std::runtime_error(exception_msg(err.str()));
  //      }
  //      for(int q = bonds[b].offsets.size() - 1; q >= 0; q--){
  //        if(idxs[b] < bonds[b].offsets[q])
  //          continue;
  //        Qidxs[b] = q;
  //        break;
  //      }
  //    }
  //    std::vector<int> Q_acc(bondNum, 1);
  //    for(int b = bondNum	- 1; b > 0; b--)
  //      Q_acc[b - 1] = Q_acc[b] * bonds[b].Qnums.size();
  //    int Qoff = 0;
  //    for(int b = 0; b < bondNum; b++)
  //      Qoff += Q_acc[b] * Qidxs[b];

  //    if(QidxEnc.find(Qoff) != QidxEnc.end()){
  //      int Q_RQoff = Qoff / CQdim;
  //      int Q_CQoff = Qoff % CQdim;
  //      Block* blk = RQidx2Blk.find(Q_RQoff)->second;
  //      size_t B_cDim = blk->Cnum;
  //      size_t sB_cDim = CQidx2Dim.find(Q_CQoff)->second;
  //      size_t blkRoff = RQidx2Off.find(Q_RQoff)->second;
  //      size_t blkCoff = CQidx2Off.find(Q_CQoff)->second;
  //      Real* boff = blk->m_elem + (blkRoff * B_cDim) + blkCoff;
  //      int cnt = 0;
  //      std::vector<int> D_acc(bondNum, 1);
  //      for(int b = bondNum	- 1; b > 0; b--)
  //        D_acc[b - 1] = D_acc[b] * bonds[b].Qdegs[Qidxs[b]];
  //      for(int b = 0; b < bondNum; b++)
  //        cnt += (idxs[b] - bonds[b].offsets[Qidxs[b]]) * D_acc[b];
  //      return boff[(cnt / sB_cDim) * B_cDim + cnt % sB_cDim];
  //    }
  //    else{
  //      return 0.0;
  //    }
  //  }
  //  catch(const std::exception& e){
  //    propogate_exception(e, "In function UniTensor::at(std::vector<size_t>&):");
  //    return 0;
  //  }
  //}

  ///*********************  Private **********************/

  ///************* developping *************/

  //Real UniTensor::norm(rflag tp) const{
  //  try{
  //    throwTypeError(tp);
  //    return vectorNorm(elem, elemNum(), 1, ongpu);
  //  }
  //  catch(const std::exception& e){
  //    propogate_exception(e, "In function UniTensor::norm(uni10::rflag ):");
  //    return 0;
  //  }
  //}

  //Real UniTensor::max(rflag tp) const{
  //  try{
  //    throwTypeError(tp);
  //    return elemMax(elem, elemNum(), ongpu);
  //  }
  //  catch(const std::exception& e){
  //    propogate_exception(e, "In function UniTensor::max(uni10::rflag ):");
  //    return 0;
  //  }
  //}

  //Real UniTensor::absMax(rflag tp) const{
  //  try{
  //    throwTypeError(tp);
  //    return elemAbsMax(elem, elemNum(), ongpu);
  //  }
  //  catch(const std::exception& e){
  //    propogate_exception(e, "In function UniTensor::absMax(uni10::rflag ):");
  //    return 0;
  //  }
  //}

  //UniTensor& UniTensor::normalize(rflag tp){
  //  try{
  //    throwTypeError(tp);
  //    Real norm = vectorNorm(elem, elemNum(), 1, ongpu);
  //    vectorScal((1./norm), elem, elemNum(), ongpu);
  //  }
  //  catch(const std::exception& e){
  //    propogate_exception(e, "In function UiTensor::normalize(uni10::rflag ):");
  //  }
  //  return *this;
  //}

  //UniTensor& UniTensor::maxNorm(rflag tp){
  //  try{
  //    throwTypeError(tp);
  //    Real max = elemMax(elem, elemNum(), ongpu);
  //    vectorScal((1./max), elem, elemNum(), ongpu);
  //  }
  //  catch(const std::exception& e){
  //    propogate_exception(e, "In function UiTensor::maxNorm(uni10::rflag ):");
  //  }
  //  return *this;
  //}

  //UniTensor& UniTensor::absMaxNorm(rflag tp){
  //  try{
  //    throwTypeError(tp);
  //    Real absMax = elemAbsMax(elem, elemNum(), ongpu);
  //    vectorScal((1./absMax), elem, elemNum(), ongpu);
  //  }
  //  catch(const std::exception& e){
  //    propogate_exception(e, "In function UiTensor::absMaxNorm(uni10::rflag ):");
  //  }
  //  return *this;
  //}

  //std::vector<UniTensor> UniTensor::hosvd(rflag tp, int* _group_labels, int* _groups, size_t _groupsSize, std::vector<Matrix>& Ls)const{
  //  try{
  //    std::vector<int> group_labels(_group_labels, _group_labels+this->bondNum());
  //    std::vector<int> groups(_groups, _groups+_groupsSize);
  //    return hosvd(tp, group_labels, groups, Ls);
  //  }
  //  catch(const std::exception& e){
  //    propogate_exception(e, "In function UniTensor::hosvd(uni10::rflag ,int* ,int* ,size_t ,std::vector<Matrix>& Ls)const;");
  //    return std::vector<UniTensor>();
  //  }
  //}

  //std::vector<UniTensor> UniTensor::hosvd(rflag tp, std::vector<int>& group_labels, std::vector<int>& groups, std::vector<Matrix>& Ls)const{
  //  try{
  //    bool withoutSymmetry = true;
  //    for(size_t b = 0; b < bonds.size(); b++){
  //      if(bonds[b].Qnums.size() != 1)
  //        withoutSymmetry = false;
  //    }
  //    if(!withoutSymmetry){
  //      std::ostringstream err;
  //      err<<"The tensor has symmetry quantum numbers. Cannot use non-symmetry version hosvd(size_t, std::vector<uni10::Matrix>&)";
  //      err<<"\n  Hint: Use UniTensor::hosvd(uni10::rflag, std::vector<int>&, std::vector<int>& , std::vector<uni10::Matrix> >&)";
  //      throw std::runtime_error(exception_msg(err.str()));
  //    }
  //    std::vector<std::map<Qnum, Matrix> > symLs;
  //    const std::vector<UniTensor>& outs = hosvd(tp , group_labels, groups, symLs, true);
  //    Ls.clear();
  //    Qnum q0(0);
  //    for(size_t i = 0; i < symLs.size(); i++)
  //      Ls.push_back(symLs[i][q0]);
  //    return outs;
  //  }
  //  catch(const std::exception& e){
  //    propogate_exception(e, "In function UniTensor::hosvd(uni10::rflag, std::vector<int>, std::vector<int>, std::vector<Matrix>&):");
  //    return std::vector<UniTensor>();
  //  }
  //}

  //std::vector<UniTensor> UniTensor::hosvd(rflag tp, int* _group_labels, int* _groups, size_t _groupsSize, std::vector<std::map<Qnum, Matrix> >& Ls, bool returnL)const{
  //  try{
  //    std::vector<int> group_labels(_group_labels, _group_labels+this->bondNum());
  //    std::vector<int> groups(_groups, _groups+_groupsSize);
  //    return hosvd(tp, group_labels, groups, Ls, returnL);
  //  }
  //  catch(const std::exception& e){
  //    propogate_exception(e, "In function UniTensor::hosvd(uni10::rflag ,int* ,int* ,size_t ,std::vector<std::map<Qnum, Matrix> >& Ls, bool returnL)const;");
  //    return std::vector<UniTensor>();
  //  }
  //}

  //std::vector<UniTensor> UniTensor::hosvd(rflag tp, std::vector<int>& group_labels, std::vector<int>& groups, std::vector<std::map<Qnum, Matrix> >& Ls, bool returnL)const{
  //  throwTypeError(tp);
  //  try{
  //    if((status & HAVEBOND) == 0){
  //      std::ostringstream err;
  //      err<<"Cannot perform higher order SVD on a tensor without bonds(scalar).";
  //      throw std::runtime_error(exception_msg(err.str()));
  //    }
  //    if((status & HAVEELEM) == 0){
  //      std::ostringstream err;
  //      err<<"Cannot perform higher order SVD on a tensor before setting its elements.";
  //      throw std::runtime_error(exception_msg(err.str()));
  //    }
  //    if(group_labels.size() != this->bondNum()){
  //      std::ostringstream err;
  //      err<<"The size of Group labels is not match to the number of tensor's bond.";
  //      throw std::runtime_error(exception_msg(err.str()));
  //    }

  //    UniTensor T(*this);
  //    size_t groupElemNum=0;
  //    for(size_t n = 0; n < groups.size(); n++)
  //      groupElemNum+=groups[n];

  //    if(returnL)
  //      Ls.assign(groups.size(), std::map<Qnum, Matrix>());
  //    std::vector<UniTensor> Us;
  //    UniTensor S(T);

  //    std::vector<int> lrsp_labels = group_labels;
  //    std::vector<int> rsp_labels = group_labels;

  //    int min = *std::min_element(rsp_labels.begin(), rsp_labels.end());

  //    for(size_t m = 0; m < groups.size(); m++){
  //      int pos=0;
  //      for(size_t l = 0; l < groupElemNum; l++){
  //        if(l >= groupElemNum-groups[m])
  //          rsp_labels[pos] = lrsp_labels[l-(groupElemNum-groups[m])];
  //        else
  //          rsp_labels[pos] = lrsp_labels[l+groups[m]];
  //        pos++;
  //      }
  //      T.permute(RTYPE, lrsp_labels, groups[m]);
  //      std::vector<Bond> bonds(T.bonds.begin(), T.bonds.begin() + groups[m]);
  //      bonds.push_back(combine(bonds).dummy_change(BD_OUT));
  //      Us.push_back(UniTensor(RTYPE, bonds));
  //      for(std::map<Qnum, Block>::iterator it = T.blocks.begin(); it != T.blocks.end(); it++){
  //        std::vector<Matrix> svd = it->second.svd(RTYPE);
  //        Us[m].putBlock(it->first, svd[0]);
  //        if(returnL)
  //          Ls[m][it->first] = svd[1];
  //      }
  //      for(int c = 0; c < groups[m]; c++)
  //        Us[m].labels[c] = lrsp_labels[c];
  //      Us[m].labels[groups[m]] = min -m - 1;
  //      UniTensor UT = Us[m];
  //      S *= UT.transpose(RTYPE);
  //      lrsp_labels = rsp_labels;
  //    } 
  //    Us.push_back(S);
  //    return Us;
  //  }
  //  catch(const std::exception& e){
  //    propogate_exception(e, "In function UniTensor::hosvd(uni10::rflag , std::vector<int>& , std::vector<int>& , std::vector<std::map<Qnum, Matrix> >& , bool returnL)const:");
  //    return std::vector<UniTensor>();
  //  }
  //}

  //std::vector<UniTensor> UniTensor::hosvd(rflag tp, size_t modeNum, size_t fixedNum)const{
  //  try{
  //    std::vector<std::map<Qnum, Matrix> > symLs;
  //    std::vector<int> group_labels=this->labels;
  //    std::vector<int> groups(modeNum, (this->bondNum()-fixedNum)/modeNum);
  //    return hosvd(tp, group_labels, groups, symLs, false);
  //  }
  //  catch(const std::exception& e){
  //    propogate_exception(e, "In function UniTensor::hosvd(uni10::rflag, size_t, size_t = 0):");
  //    return std::vector<UniTensor>();
  //  }
  //}

  //std::vector<UniTensor> UniTensor::hosvd(rflag tp, size_t modeNum, size_t fixedNum, std::vector<std::map<Qnum, Matrix> >& Ls)const{
  //  try{
  //    std::vector<std::map<Qnum, Matrix> > symLs;
  //    std::vector<int> group_labels=this->labels;
  //    std::vector<int> groups(modeNum, (this->bondNum()-fixedNum)/modeNum );
  //    return hosvd(tp, group_labels, groups, Ls, true);
  //  }
  //  catch(const std::exception& e){
  //    propogate_exception(e, "In function UniTensor::hosvd(uni10::rflag, size_t, size_t, std::vector<std::map<uni10::Qnum, uni10::Matrix> >&):");
  //    return std::vector<UniTensor>();
  //  }
  //}

  //std::vector<UniTensor> UniTensor::hosvd(rflag tp, size_t modeNum, size_t fixedNum, std::vector<Matrix>& Ls)const{
  //  try{
  //    bool withoutSymmetry = true;
  //    for(size_t b = 0; b < bonds.size(); b++){
  //      if(bonds[b].Qnums.size() != 1)
  //        withoutSymmetry = false;
  //    }
  //    if(!withoutSymmetry){
  //      std::ostringstream err;
  //      err<<"The tensor has symmetry quantum numbers. Cannot use non-symmetry version hosvd(size_t, std::vector<uni10::Matrix>&)";
  //      err<<"\n  Hint: Use UniTensor::hosvd(uni10::rflag, size_t, size_t, std::vector<std::map<uni10::Qnum, uni10::Matrix> >&)";
  //      throw std::runtime_error(exception_msg(err.str()));
  //    }
  //    std::vector<std::map<Qnum, Matrix> > symLs;
  //    const std::vector<UniTensor>& outs = hosvd(tp , modeNum, fixedNum, symLs);
  //    Ls.clear();
  //    Qnum q0(0);
  //    for(size_t i = 0; i < symLs.size(); i++)
  //      Ls.push_back(symLs[i][q0]);
  //    return outs;
  //  }
  //  catch(const std::exception& e){
  //    propogate_exception(e, "In function UniTensor::hosvd(uni10::rflag, size_t, size_t, std::vector<Matrix>&):");
  //    return std::vector<UniTensor>();
  //  }
  //}

  template class UniTensor<uni10_double64>;
  template class UniTensor<uni10_complex128>;
}

