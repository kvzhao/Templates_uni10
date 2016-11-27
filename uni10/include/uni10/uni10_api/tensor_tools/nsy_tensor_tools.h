#ifndef __UNI10_NSY_TENSOR_TOOLS_H__
#define __UNI10_NSY_TENSOR_TOOLS_H__

#include <stdio.h>

#include "uni10/uni10_api/rng.h"
#include "uni10/uni10_api/UniTensor.h"

namespace uni10{
  
  namespace tensor_tools{

    //Function prototype.
    template<typename uni10_type>
      U_para<uni10_type>* init_para_nsy(U_para<uni10_type>* para);

    template<typename uni10_type>
      void copy_para_nsy(U_para<uni10_type>* para, const U_para<uni10_type>* src_para);

    template<typename uni10_type>
      void free_para_nsy(U_para<uni10_type>* para);

    template<typename uni10_type>
      void init_nsy(U_para<uni10_type>* para);

    template<typename uni10_type>
      uni10_uint64 grouping_nsy(U_para<uni10_type>* para);

    template <typename uni10_type>
      void initBlocks_nsy(U_para<uni10_type>* para);

    template <typename uni10_type>
      void putBlock_nsy(U_para<uni10_type>* para,const Qnum& qnum, const Block<uni10_type>& mat);

    template <typename uni10_type>
      void set_zero_nsy(U_para<uni10_type>* para);

    template <typename uni10_type>
      void randomize_nsy(U_para<uni10_type>* para);

    template <typename uni10_type>
      void permute_nsy(const U_para<uni10_type>* T1_para, const std::vector<uni10_int32>& rsp_outin,
        U_para<uni10_type>* T2_para, uni10_bool inorder);

    //Functions.
    template<typename uni10_type>
       U_para<uni10_type>* init_para_nsy(U_para<uni10_type>* para){

        //std::cout << para-> check_status;
        para = new struct U_para<uni10_type>[1];
        para->check_status = 1;
        //std::cout << para-> check_status;
        //exit(0);
        para->nsy = new struct no_sym_para<uni10_type>[1];

        return para;
        
      }

    template<typename uni10_type>
      void copy_para_nsy(U_para<uni10_type>* para, const U_para<uni10_type>* src_para){

        *para->nsy = *src_para->nsy;

      }

    template<typename uni10_type>
       void free_para_nsy(U_para<uni10_type>* para){
        
         delete [] para->nsy;
         delete [] para;
        
      }

    template<typename uni10_type>
      void init_nsy(U_para<uni10_type>* para){

        if(para->nsy->bonds.size()){
          para->nsy->U_elemNum = grouping_nsy(para);
          if(!(para->nsy->blocks.size() > 0)){ //No block in Tensor, Error!
            uni10_error_msg(true, "%s", "There is no symmetry block with the given bonds:\n");
            for(uni10_int32 b = 0; b < (uni10_int32)para->nsy->bonds.size(); b++)
              std::cout<<"    "<<para->nsy->bonds[b];
          }

          para->nsy->labels.assign(para->nsy->bonds.size(), 0);
          for(uni10_int32 b = 0; b < (uni10_int32)para->nsy->bonds.size(); b++)
            para->nsy->labels[b] = b;
          para->nsy->status |= UniTensor<uni10_type>::GET_HAVEBOND();

        }
        else{
          Qnum q0(0);
          para->nsy->blocks[q0] = Block<uni10_type>(1, 1);
          para->nsy->RBondNum = 0;
          para->nsy->Rdim = 0;
          para->nsy->Cdim = 0;
          para->nsy->U_elemNum = 1;
          para->nsy->status |= UniTensor<uni10_type>::GET_HAVEELEM();
        }

        //UniTensor<uni10_type>::GET_ELEMNUM()+= para->nsy->U_elemNum;
        //UniTensor<uni10_type>::GET_COUNTER()++;
        //if((uni10_uint64)UniTensor<uni10_type>::ELEMNUM > UniTensor<uni10_type>::GET_MAXELEMNUM())
        //  UniTensor<uni10_type>::GET_MAXELEMNUM() = UniTensor<uni10_type>::GET_ELEMNUM();
        //if(para->nsy->U_elemNum > UniTensor<uni10_type>::GET_MAXELEMTEN())
        //  UniTensor<uni10_type>::GET_MAXELEMTEN() = para->nsy->U_elemNum;

        para->nsy->U_elem.init(1, para->nsy->U_elemNum, false);
        initBlocks_nsy(para);

      }

    template<typename uni10_type>
      uni10_uint64 grouping_nsy(U_para<uni10_type>* para){

        para->nsy->blocks.clear();
        Qnum q0(0);
        uni10_int32   row_bondNum = 0;
        uni10_int32   col_bondNum = 0;
        para->nsy->Rdim = 1;
        para->nsy->Cdim = 1;
        bool IN_BONDS_BEFORE_OUT_BONDS = true;
        for(uni10_int32 i = 0; i < (uni10_int32)para->nsy->bonds.size(); i++){
          uni10_error_msg( para->nsy->bonds[i].const_getQdegs().size() > 1 || para->nsy->bonds[i].const_getQnums().size() > 1, "%s", "Ihis UniTensor has symmetry!!");
          if(para->nsy->bonds[i].type() == BD_IN){
            uni10_error_msg(!(IN_BONDS_BEFORE_OUT_BONDS == true), 
                "%s","Error in the input bond array: BD_OUT bonds must be placed after all BD_IN bonds.");
            para->nsy->Rdim *= para->nsy->bonds[i].const_getQdegs()[0];
            row_bondNum++;
          }
          else{
            para->nsy->Cdim *= para->nsy->bonds[i].const_getQdegs()[0];
            col_bondNum++;
            IN_BONDS_BEFORE_OUT_BONDS = false;
          }
        }
        para->nsy->RBondNum = row_bondNum;
        para->nsy->blocks[q0] = Block<uni10_type>(para->nsy->Rdim, para->nsy->Cdim);

        return para->nsy->Rdim*para->nsy->Cdim;

      }

    template <typename uni10_type>
      void initBlocks_nsy(U_para<uni10_type>* para){
        uni10_uint64 offset = 0;
        typename std::map< Qnum, Block<uni10_type> >::iterator it = para->nsy->blocks.begin();
        for(; it != para->nsy->blocks.end(); it++ ){
          it->second.elem_enforce().__elem = &(para->nsy->U_elem.__elem[offset]);
          offset += it->second.row_enforce() * it->second.col_enforce();
        }
      }

    template <typename uni10_type>
      void putBlock_nsy(U_para<uni10_type>* para,const Qnum& qnum, const Block<uni10_type>& mat){

        typename std::map<Qnum, Block<uni10_type> >::iterator it;

        if(!((it = para->nsy->blocks.find(qnum)) != para->nsy->blocks.end())){
          uni10_error_msg(true, "%s", "There is no block with the given quantum number ");
          std::cout<<qnum;
        }

        uni10_error_msg(!(mat.row() == it->second.row_enforce() && mat.col() == it->second.col_enforce()), "%s", 
            "The dimension of input matrix does not match for the dimension of the block with quantum number \n  Hint: Use Matrix::resize(int, int)");

        if(mat.getElem() != it->second.getElem()){
          if(mat.isDiag()){
            uni10_error_msg(true, "%s","Developping!!!");
          }
          else
            it->second.elem_enforce().copy(0, mat.const_elem_enforce(), it->second.row_enforce() * it->second.col_enforce() );
        }

        para->nsy->status |= UniTensor<uni10_type>::GET_HAVEELEM();

      }

    template <typename uni10_type>
      void set_zero_nsy(U_para<uni10_type>* para){

        para->nsy->U_elem.set_zeros();

      }

    template <typename uni10_type>
      void randomize_nsy(U_para<uni10_type>* para){

        uni10_error_msg(true, "%s", "Developping!!!\n");

      }

    template <typename uni10_type>
      void permute_nsy(const U_para<uni10_type>* T1_para, const std::vector<uni10_int32>& rsp_outin,
          U_para<uni10_type>* T2_para, uni10_bool inorder){

        uni10_uint64 bondNum = T1_para->nsy->bonds.size();
        if(!inorder){
          std::vector<size_t> transAcc(bondNum);
          std::vector<size_t> newAcc(bondNum);
          transAcc[bondNum - 1] = 1;
          newAcc[bondNum - 1] = 1;
          for(int b = bondNum - 1; b > 0; b--)
            newAcc[b - 1] = newAcc[b] * T2_para->nsy->bonds[b].const_getQdegs()[0];
          std::vector<int> bondDims(bondNum);
          std::vector<int> idxs(bondNum);
          for(int b = 0; b < (int)bondNum; b++){
            transAcc[rsp_outin[b]] = newAcc[b];
            bondDims[b] = T1_para->nsy->bonds[b].const_getQdegs()[0];
            idxs[b] = 0;
          }
          size_t cnt_ot = 0;
          for(size_t i = 0; i < T2_para->nsy->U_elemNum; i++){
            T2_para->nsy->U_elem.__elem[cnt_ot] = T1_para->nsy->U_elem.__elem[i];
            for(int bend = bondNum - 1; bend >= 0; bend--){
              idxs[bend]++;
              if(idxs[bend] < bondDims[bend]){
                cnt_ot += transAcc[bend];
                break;
              }
              else{
                cnt_ot -= transAcc[bend] * (idxs[bend] - 1);
                idxs[bend] = 0;
              }
            }
          }
        }
        else{  //non-symmetry inorder
          T2_para->nsy->U_elem.copy(0, T1_para->nsy->U_elem, T1_para->nsy->U_elemNum);
        }
        
      }

  };

};

#endif
