#ifndef __UNI10_NSY_TENSOR_TOOLS_H__
#define __UNI10_NSY_TENSOR_TOOLS_H__

#include <stdio.h>

#include "uni10/uni10_api/UniTensor.h"

namespace uni10{
  
  namespace tensor_tools{

    //Function prototype.
    template<typename uni10_type>
      U_para<uni10_type>* init_para_nsy(U_para<uni10_type>* para);

    template<typename uni10_type>
       void free_para_nsy(U_para<uni10_type>* para);

    template<typename uni10_type>
      void init_nsy(U_para<uni10_type>* para);

    template<typename uni10_type>
      uni10_uint64 grouping_nsy(U_para<uni10_type>* para);

    template <typename uni10_type>
      void initBlocks_nsy(U_para<uni10_type>* para);

    //Functions.
    template<typename uni10_type>
       U_para<uni10_type>* init_para_nsy(U_para<uni10_type>* para){

        para = new struct U_para<uni10_type>[1];
        para->nsy = new struct no_sym_para<uni10_type>[1];

        return para;
        
      }

    template<typename uni10_type>
       void free_para_nsy(U_para<uni10_type>* para){
        
         delete para->nsy;
         delete para;
        
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

  };

};

#endif
