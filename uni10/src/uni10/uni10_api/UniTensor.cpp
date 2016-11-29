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
#include "uni10/uni10_api/linalg.h"
#include "uni10/uni10_api/tensor_tools/tensor_tools.h"

namespace uni10{

  template <typename uni10_type>
    UniTensor<uni10_type>::UniTensor(): style(no_sym), paras(NULL){

    };

  template <typename uni10_type>
    UniTensor<uni10_type>::UniTensor(uni10_type val){ 

      style = check_bonds();
      this->init_para();
      this->meta_link();
      *status = 0;
      this->init();
      this->U_elem->__elem[0] = val;
    }

  template <typename uni10_type>
    UniTensor<uni10_type>::UniTensor(const std::vector<Bond>& _bonds, const std::string& _name){
     
      style = check_bonds();
      this->init_para();
      this->meta_link();
      *name  = _name;
      *bonds = _bonds;
      *status= 0;
      this->init();
       
    }

  template <typename uni10_type>
    UniTensor<uni10_type>::UniTensor(const std::vector<Bond>& _bonds, int* _labels, const std::string& _name){

      style = check_bonds();
      this->init_para();
      this->meta_link();
      *name  = _name;
      *bonds = _bonds;
      *status= 0;
      this->init();
      this->setLabel(_labels);

    }

  template <typename uni10_type>
    UniTensor<uni10_type>::UniTensor(const std::vector<Bond>& _bonds, std::vector<int>& _labels, const std::string& _name){

      style = check_bonds();
      this->init_para();
      this->meta_link();
      *name  = _name;
      *bonds = _bonds;
      *status= 0;
      this->init();
      this->setLabel(_labels);

    }

  template <typename uni10_type>
    UniTensor<uni10_type>::UniTensor(const UniTensor& UniT): style(UniT.style){

      this->init_para();
      this->meta_link();
      this->copy_para(UniT.paras);
      this->initBlocks();
     
    }

  //template <typename uni10_type>
  //  UniTensor<uni10_type>::UniTensor(const std::string& fname){

  //  }

  template <typename uni10_type>
    UniTensor<uni10_type>::UniTensor(const Block<uni10_type>& blk){

      this->style = no_sym;
      Bond bdi(BD_IN, blk.Rnum);
      Bond bdo(BD_OUT, blk.Cnum);
      this->init_para();
      this->meta_link();
      bonds->push_back(bdi);
      bonds->push_back(bdo);
      this->init();
      this->setElem(blk.getElem());
      uni10_error_msg(true, "%s", "Developping!!\n");
      //this->putBlock(blk);

    }

  template <typename uni10_type>
    UniTensor<uni10_type>::~UniTensor(){
      UniTensor<uni10_type>::ELEMNUM -= *U_elemNum;
      UniTensor<uni10_type>::COUNTER--;
      free_para();
    }

  template <typename uni10_type>
    UniTensor<uni10_type>& UniTensor<uni10_type>::assign(const std::vector<Bond>& _bond ){
  
      UniTensor<uni10_type> T(_bond);
      *this = T;
      return *this;

    }

  template <typename uni10_type>
    void UniTensor<uni10_type>::putBlock(const Block<uni10_type>& mat){
      Qnum q0(0);
      putBlock(q0, mat);
    }

  template <typename uni10_type>
    void UniTensor<uni10_type>::putBlock(const Qnum& qnum, const Block<uni10_type>& mat){

      tensor_tools::putBlock(paras, qnum, mat, style);

    }

  template <typename uni10_type>
    void UniTensor<uni10_type>::setName(const std::string& _name){
      *name = _name;
    }

  template <typename uni10_type>
    void UniTensor<uni10_type>::setLabel(const uni10_int32 newLabel, const uni10_uint64 idx){
      uni10_error_msg(labels->size() <= idx, "%s", "The bond index is out of the range of vector(labels).");
      (*labels)[idx] = newLabel;
    }

  template <typename uni10_type>
    void UniTensor<uni10_type>::setLabel(const std::vector<uni10_int32>& newLabels){
      uni10_error_msg(!(bonds->size() == newLabels.size()), "%s", "The size of input vector(labels) does not match for the number of bonds.");
      *labels = newLabels;
    }

  template <typename uni10_type>
    void UniTensor<uni10_type>::setLabel(uni10_int32* newLabels){
      std::vector<uni10_int32> labels(newLabels, newLabels + bonds->size());
      setLabel(labels);
    }

  template <typename uni10_type>
    const Block<uni10_type>& UniTensor<uni10_type>::const_getBlock()const{
      Qnum q0(0);
      return const_getBlock(q0);
    }

  template <typename uni10_type>
    const Block<uni10_type>& UniTensor<uni10_type>::const_getBlock(const Qnum& qnum)const{
      typename std::map<Qnum, Block<uni10_type> >::const_iterator it = blocks->find(qnum);
      if(it == blocks->end()){
        uni10_error_msg(true, "%s", "There is no block with the given quantum number ");
        std::cout << qnum;
      }
      return it->second;
    }

  template <typename uni10_type>
    std::vector<Qnum> UniTensor<uni10_type>::blockQnum()const{
      std::vector<Qnum> keys;
      typename std::map<Qnum, Block<uni10_type> >::const_iterator it = blocks->begin();
      for(; it != blocks->end(); it++)
        keys.push_back(it->first);
      return keys;
    }

  template <typename uni10_type>
    Qnum UniTensor<uni10_type>::blockQnum(uni10_uint64 idx)const{

      uni10_error_msg(!(idx < blocks->size()), "Index exceeds the number of the blocks( %ld ).", blocks->size());
      typename std::map<Qnum, Block<uni10_type> >::const_iterator it = blocks->begin();
      for(; it != blocks->end(); it++){
        if(idx == 0)
          return it->first;
        idx--;
      }
      return Qnum();
    }

  template <typename uni10_type>
    std::map< Qnum, Matrix<uni10_type> > UniTensor<uni10_type>::getBlocks()const{
      std::map<Qnum, Matrix<uni10_type> > mats;
      typename std::map<Qnum, Block<uni10_type> >::const_iterator it = blocks->begin();
      for(; it != blocks->end(); it++){
        Matrix<uni10_type> mat(it->second.Rnum, it->second.Cnum, it->second.elem.__elem);
        mats.insert(std::pair<Qnum, Matrix<uni10_type> >(it->first, mat));
      }
      return std::map< Qnum, Matrix<uni10_type> >();
    }

  template <typename uni10_type>
    Matrix<uni10_type> UniTensor<uni10_type>::getBlock(uni10_bool diag)const{
      Qnum q0(0);
      return getBlock(q0, diag);
    }

  template <typename uni10_type>
    Matrix<uni10_type> UniTensor<uni10_type>::getBlock(const Qnum& qnum, uni10_bool diag)const{
      typename std::map<Qnum, Block<uni10_type> >::const_iterator it = blocks->find(qnum);
      if(it == blocks->end()){
        uni10_error_msg(true, "%s", "There is no block with the given quantum number ");
        std::cout<<qnum;
      }
      if(diag)
        return getDiag(it->second);
      else{
        Matrix<uni10_type> mat(it->second.Rnum, it->second.Cnum, it->second.elem.__elem, false);
        return mat;
      }
    }

  template <typename uni10_type>
    void UniTensor<uni10_type>::setElem(const uni10_type* _elem){
      UELEM(uni10_elem, _package, _type)<uni10_type> _src(_elem, 1, *U_elemNum, false);
      U_elem->copy(0, _src, *U_elemNum);
      *status |= HAVEELEM;
  }

  template <typename uni10_type>
    void UniTensor<uni10_type>::setElem(const std::vector<uni10_type>& _elem){
      uni10_error_msg(_elem.size() != *U_elemNum, "%s", "The number of input elements is defferent from the size of the UniTensor");
      setElem(&_elem[0]);
    }

  template <typename uni10_type>
    void UniTensor<uni10_type>::addGate(const std::vector<_Swap>& swaps){
      uni10_error_msg((*status & HAVEBOND) == 0, "%s", "Adding swap gates to a tensor without bonds(scalar).");
      uni10_error_msg((*status & HAVEELEM) == 0,"%s", "Cannot add swap gates to a tensor before setting its elements.")

        int sign = 1;
      int bondNum = bonds->size();
      std::vector<int> Q_idxs(bondNum, 0);
      std::vector<int> Q_Bdims(bondNum, 0);
      for(int b = 0; b < bondNum; b++)
        Q_Bdims[b] = bonds[b].Qnums.size();
      int Q_off;
      int tmp;
      int RQoff, CQoff;
      size_t sB_r, sB_c;	//sub-block of a Qidx
      size_t sB_rDim, sB_cDim;	//sub-block of a Qidx
      size_t B_cDim;
      uni10_type* Eptr;
      typename std::map<int, size_t>::iterator it = QidxEnc.begin()
    for(; it != QidxEnc.end(); it++){
      Q_off = it->first;
      tmp = Q_off;
      for(int b = bondNum - 1; b >= 0; b--){
        Q_idxs[b] = tmp % Q_Bdims[b];
        tmp /= Q_Bdims[b];
      }
      RQoff = Q_off / CQdim;
      CQoff = Q_off % CQdim;
      B_cDim = RQidx2Blk[RQoff]->Cnum;
      Eptr = RQidx2Blk[RQoff]->m_elem + (RQidx2Off[RQoff] * B_cDim) + CQidx2Off[CQoff];
      sB_rDim = RQidx2Dim[RQoff];
      sB_cDim = CQidx2Dim[CQoff];

      int sign01 = 0;
      for(size_t i = 0; i < swaps.size(); i++)
        sign01 ^= (bonds[swaps[i].b1].Qnums[Q_idxs[swaps[i].b1]].prtF() & bonds[swaps[i].b2].Qnums[Q_idxs[swaps[i].b2]].prtF());
      sign = sign01 ? -1 : 1;

      for(sB_r = 0; sB_r < sB_rDim; sB_r++)
        for(sB_c = 0; sB_c < sB_cDim; sB_c++)
          Eptr[(sB_r * B_cDim) + sB_c] *= sign;
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function UniTensor::addGate(uni10::rflag, std::vector<_Swap>&):");
  }

    }

  template <typename uni10_type>
    std::vector<_Swap> exSwap(const UniTensor<uni10_type>& Tb)const{

    }


  template <typename uni10_type>
    void UniTensor<uni10_type>::printDiagram()const{
      if(!(*status & HAVEBOND)){
        if(U_elem->__ongpu)
          std::cout<<"\nScalar: " << U_elem->__elem[0]<<", onGPU";
        else
          std::cout<<"\nScalar: " << U_elem->__elem[0];
        std::cout<<"\n\n";
      }
      else{

        uni10_uint64 row = 0;
        uni10_uint64 col = 0;

        std::vector<Bond>_bonds = *bonds;
        for(uni10_uint64 i = 0; i < _bonds.size(); i++)
          if(_bonds[i].type() == BD_IN)
            row++;
          else
            col++;
        uni10_uint64 layer = std::max(row, col);
        uni10_uint64 nmlen = (*name).length() + 2;
        uni10_uint64 star = 12 + (14 - nmlen) / 2;
        for(uni10_uint64 s = 0; s < star; s++)
          std::cout << "*";
        if((*name).length() > 0)
          std::cout << " " << *name << " ";
        for(uni10_uint64 s = 0; s < star; s++)
          std::cout <<"*";
        std::cout<<std::endl;

        if(U_elem->__uni10_typeid == 1)
          std::cout << "REAL" << std::endl;
        else if(U_elem->__uni10_typeid == 2)
          std::cout << "COMPLEX" << std::endl;

        if(U_elem->__ongpu)
          std::cout<<"\n                 onGPU";
        std::cout << "\n             ____________\n";
        std::cout << "            |            |\n";
        uni10_uint64 llab = 0;
        uni10_uint64 rlab = 0;
        char buf[128];
        for(uni10_uint64 l = 0; l < layer; l++){
          if(l < row && l < col){
            llab = (*labels)[l];
            rlab = (*labels)[row + l];
            sprintf(buf, "    %5ld___|%-4d    %4d|___%-5ld\n", llab, _bonds[l].dim(), _bonds[row + l].dim(), rlab);
            std::cout<<buf;
          }
          else if(l < row){
            llab = (*labels)[l];
            sprintf(buf, "    %5ld___|%-4d    %4s|\n", llab, _bonds[l].dim(), "");
            std::cout<<buf;
          }
          else if(l < col){
            rlab = (*labels)[row + l];
            sprintf(buf, "    %5s   |%4s    %4d|___%-5ld\n", "", "", _bonds[row + l].dim(), rlab);
            std::cout << buf;
          }
          std::cout << "            |            |   \n";
        }
        std::cout << "            |____________|\n";

        std::cout << "\n================BONDS===============\n";
        for(uni10_uint64 b = 0; b < _bonds.size(); b++)
          std::cout << _bonds[b];

        std::cout << "\n\nTotal elemNum: "<<(*U_elemNum)<<std::endl;
        std::cout << "====================================\n";
      }
    }

  template <typename uni10_type>
    void UniTensor<uni10_type>::init_para(){

      paras = tensor_tools::init_para(paras, style);

    }

  template <typename uni10_type>
    void UniTensor<uni10_type>::copy_para(U_para<uni10_type>* src_para){

      tensor_tools::copy_para(paras, src_para, style);

    }

  template <typename uni10_type>
    void UniTensor<uni10_type>::meta_link(){

      if(style == 0){

        name      = &paras->nsy->name;
        bonds     = &paras->nsy->bonds;
        labels    = &paras->nsy->labels;
        RBondNum  = &paras->nsy->RBondNum;  
        Rdim      = &paras->nsy->Rdim;
        Cdim      = &paras->nsy->Cdim;
        U_elemNum = &paras->nsy->U_elemNum;    
        blocks    = &paras->nsy->blocks;
        U_elem    = &paras->nsy->U_elem;
        status    = &paras->nsy->status;

      }
      else if(style == 1){
        name      = &paras->bsy->name;
        bonds     = &paras->bsy->bonds;
        labels    = &paras->bsy->labels;
        RBondNum  = &paras->bsy->RBondNum;  
        blocks    = &paras->bsy->blocks;
        uni10_error_msg(true, "%s", "Developping!!!");

      }
      else if(style == 2){
        //name = &paras->ssy->name;
        uni10_error_msg(true, "%s", "Developping!!!");
      }

    }

  template <typename uni10_type>
    void UniTensor<uni10_type>::init(){
      // You should init_para() first. Then you can use this function to initialize the UniTensor.
      tensor_tools::init(paras, style);
    }

  template <typename uni10_type>
    void UniTensor<uni10_type>::initBlocks(){

      tensor_tools::initBlocks(paras, style);

    }


  template <typename uni10_type>
    void UniTensor<uni10_type>::free_para(){
      // You should init_para() first. Then you can use this function to initialize the UniTensor.
      tensor_tools::free_para(paras, style);
    }

  template <typename uni10_type> 
    contain_type UniTensor<uni10_type>::check_bonds()const{
      // blk_sym && spar_sym are developping
      contain_type s = no_sym;
      return s;
    }

  template class UniTensor<uni10_double64>;
  template class UniTensor<uni10_complex128>;
}

