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
      void copy_para_bsy(U_para<uni10_type>* para, const U_para<uni10_type>* src_para);

    template<typename uni10_type>
       void free_para_bsy(U_para<uni10_type>* para);

    template<typename uni10_type>
      void init_bsy(U_para<uni10_type>* para);

    template<typename uni10_type>
      uni10_uint64 grouping_bsy(U_para<uni10_type>* _para);

    template <typename uni10_type>
      void initBlocks_bsy(U_para<uni10_type>* para);

    template <typename uni10_type>
      void putBlock_bsy(U_para<uni10_type>* para,const Qnum& qnum, const Block<uni10_type>& mat);

    template <typename uni10_type>
      void set_zero_bsy(U_para<uni10_type>* para);

    template <typename uni10_type>
      void randomize_bsy(U_para<uni10_type>* para);

    template <typename uni10_type>
      void permute_bsy(const U_para<uni10_type>* T1_para, const std::vector<uni10_int32>& rsp_outin,
        U_para<uni10_type>* T2_para, uni10_bool inorder);

    // Functions.
    template<typename uni10_type>
      U_para<uni10_type>* init_para_bsy(U_para<uni10_type>* para){

        para = new struct U_para<uni10_type>[1];
        para->bsy = new struct blk_sym_para<uni10_type>[1];
      
        return para;   
      }

    template<typename uni10_type>
      void copy_para_bsy(U_para<uni10_type>* para, const U_para<uni10_type>* src_para){

        *para->bsy = *src_para->bsy;

      }

    template<typename uni10_type>
       void free_para_bsy(U_para<uni10_type>* para){
        
         delete para->bsy;
         delete para;
        
      }

    template<typename uni10_type>
      void init_bsy(U_para<uni10_type>* _para){

        uni10_error_msg(true, "%s", "Developping");
      }

    template<typename uni10_type>
      uni10_uint64 grouping_bsy(U_para<uni10_type>* _para){

        uni10_error_msg(true, "%s", "Developping");

      }

    template <typename uni10_type>
      void initBlocks_bsy(U_para<uni10_type>* para){

        uni10_error_msg(true, "%s", "Developping");

      }

    template <typename uni10_type>
      void putBlock_bsy(U_para<uni10_type>* para,const Qnum& qnum, const Block<uni10_type>& mat){

        uni10_error_msg(true, "%s", "Developping");

      }

    template <typename uni10_type>
      void set_zero_bsy(U_para<uni10_type>* para){

        uni10_error_msg(true, "%s", "Developping");

      }

    template <typename uni10_type>
      void randomize_bsy(U_para<uni10_type>* para){

        uni10_error_msg(true, "%s", "Developping");

      }

    template <typename uni10_type>
      void permute_bsy(const U_para<uni10_type>* T1_para, const std::vector<uni10_int32>& rsp_outin,
          U_para<uni10_type>* T2_para, uni10_bool inorder){

        uni10_error_msg(true, "%s", "Developping");

      }

  };

};

#endif
