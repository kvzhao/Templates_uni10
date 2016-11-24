#ifndef __UNI10_NSY_TENSOR_TOOLS_H__
#define __UNI10_NSY_TENSOR_TOOLS_H__

#include <stdio.h>

#include "uni10/uni10_api/UniTensor.h"

namespace uni10{
  
  namespace tensor_tools{

    //Function prototype.
    template<typename uni10_type>
      void init_nsy(U_para<uni10_type>* para);

    template<typename uni10_type>
      void grouping_nsy(U_para<uni10_type>* _para);

    template <typename uni10_type>
      void initBlocks_nsy(U_para<uni10_type>* para);

    template<typename uni10_type>
      void init_val_nsy(uni10_type val, U_para<uni10_type>* para);

    //Functions.
    template<typename uni10_type>
      void init_nsy(U_para<uni10_type>* para){

        para->nsy = new struct nsy_sym_para[1];
        para->nsy->status = 0;

        if(para->bsy->bonds.size()){
          para->bsy->U_elemNum = grouping_nsy(para);
          if(!(para->bsy->blocks.size() > 0)){ //No block in Tensor, Error!
            uni10_error_msg(true, "%s", "There is no symmetry block with the given bonds:\n");
            for(uni10_int32 b = 0; b < (uni10_int32)para->bsy->bonds.size(); b++)
              std::cout<<"    "<<para->bsy->bonds[b];
          }

          para->bsy->labels.assign(para->bsy->bonds.size(), 0);
          for(uni10_int32 b = 0; b < (uni10_int32)para->bsy->bonds.size(); b++)
            para->bsy->labels[b] = b;
          para->bsy->status |= UniTensor<uni10_type>::HAVEBOND;

        }
        else{
          Qnum q0(0);
          para->bsy->blocks[q0] = Block<uni10_type>(1, 1);
          para->bsy->RBondNum = 0;
          para->bsy->Rdim = 0;
          para->bsy->Cdim = 0;
          para->bsy->U_elemNum = 1;
          para->bsy->status |= UniTensor<uni10_type>::HAVEELEM;
        }

        UniTensor<uni10_type>::ELEMNUM += para->bsy->U_elemNum;
        UniTensor<uni10_type>::COUNTER++;
        if((uni10_uint64)UniTensor<uni10_type>::ELEMNUM > UniTensor<uni10_type>::MAXELEMNUM)
          UniTensor<uni10_type>::MAXELEMNUM = UniTensor<uni10_type>::ELEMNUM;
        if(para->bsy->U_elemNum > UniTensor<uni10_type>::MAXELEMTEN)
          UniTensor<uni10_type>::MAXELEMTEN = para->bsy->U_elemNum;

        para->bsy->U_elem.init(1, para->bsy->U_elemNum, false);
        initBlocks(para);

      }

    template<typename uni10_type>
      void grouping_nsy(U_para<uni10_type>* para){

        para->nsy->blocks.clear();
        Qnum q0(0);
        uni10_int32   row_bondNum = 0;
        uni10_int32   col_bondNum = 0;
        para->nsy->Rdim = 1;
        para->nsy->Cdim = 1;
        bool IN_BONDS_BEFORE_OUT_BONDS = true;
        for(uni10_int32 i = 0; i < (uni10_int32)para->nsy->bonds.size(); i++){
          uni10_error_msg( para->nsy->bonds[i].Qdegs.size() > 1 || para->nsy->bonds[i].Qnums.size() > 1, "%s", "Ihis UniTensor has symmetry!!");
          if(para->nsy->bonds[i].type() == BD_IN){
            uni10_error_msg(!(IN_BONDS_BEFORE_OUT_BONDS == true), 
                "%s","Error in the input bond array: BD_OUT bonds must be placed after all BD_IN bonds.");
            para->nsy->Rdim *= para->nsy->bonds[i].Qdegs[0];
            row_bondNum++;
          }
          else{
            para->nsy->Cdim *= para->nsy->bonds[i].Qdegs[0];
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
          it->second.elem.__elem = &(para->nsy->U_elem.__elem[offset]);
          offset += it->second.Rnum * it->second.Cnum;
        }
      }

    template<typename uni10_type>
      void init_val_nsy(uni10_type val, U_para<uni10_type>* para){
        init_nsy(para);
        para->nsy->U_elem.__elem[0] = val;

      }

  };

};

#endif
