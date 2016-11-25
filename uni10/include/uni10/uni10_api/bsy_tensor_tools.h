#ifndef __UNI10_BSY_TENSOR_TOOLS_H__
#define __UNI10_BSY_TENSOR_TOOLS_H__

#include <stdio.h>

#include "uni10/uni10_api/UniTensor.h"

namespace uni10{

  namespace tensor_tools{

    //Function prototype.
    template<typename uni10_type>
      U_para<uni10_type>* init_para_bsy(U_para<uni10_type>* para);

    template<typename uni10_type>
       void free_para_bsy(U_para<uni10_type>* para);

    template<typename uni10_type>
      void init_bsy(U_para<uni10_type>* para);

    template<typename uni10_type>
      uni10_uint64 grouping_bsy(U_para<uni10_type>* _para);

    template <typename uni10_type>
      void initBlocks_bsy(U_para<uni10_type>* para);

    // Functions.
    template<typename uni10_type>
      U_para<uni10_type>* init_para_bsy(U_para<uni10_type>* para){

        para = new struct U_para<uni10_type>[1];
        para->bsy = new struct blk_sym_para<uni10_type>[1];
      
        return para;   
      }

    template<typename uni10_type>
       void free_para_bsy(U_para<uni10_type>* para){
        
         delete para->bsy;
         delete para;
        
      }

    template<typename uni10_type>
      void init_bsy(U_para<uni10_type>* _para){


      }

    template<typename uni10_type>
      uni10_uint64 grouping_bsy(U_para<uni10_type>* _para){

      }

    template <typename uni10_type>
      void initBlocks_bsy(U_para<uni10_type>* para){

      }


  };

};

#endif
