#ifndef __UNI10_BSY_TENSOR_TOOLS_H__
#define __UNI10_BSY_TENSOR_TOOLS_H__

#include <stdio.h>

#include "uni10/uni10_api/UniTensor.h"

namespace uni10{

  namespace tensor_tools{

    //Function prototype.
    template<typename uni10_type>
      void init_bsy(U_para<uni10_type>* para);

    template<typename uni10_type>
      void grouping_bsy(U_para<uni10_type>* _para);

    template <typename uni10_type>
      void initBlocks_bsy(U_para<uni10_type>* para);

    template<typename uni10_type>
      void init_val_bsy(uni10_type val, U_para<uni10_type>* para);



    // Functions.
    template<typename uni10_type>
      void init_val_bsy(uni10_type val, U_para<uni10_type>* para){

      }

    template<typename uni10_type>
      void init_bsy(U_para<uni10_type>* _para){


      }

    template<typename uni10_type>
      void grouping_bsy(U_para<uni10_type>* _para){

      }

    template <typename uni10_type>
      void initBlocks_bsy(U_para<uni10_type>* para){

      }


  };

};

#endif
