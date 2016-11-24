#ifndef __UNI10_TENSOR_TOOLS_H__
#define __UNI10_TENSOR_TOOLS_H__

#include <stdio.h>

#include "uni10/uni10_api/nsy_tensor_tools.h"
#include "uni10/uni10_api/bsy_tensor_tools.h"

namespace uni10{

  namespace tensor_tools{
    // The function tools to initialize UniTensor.
    static void (*init_d[])(U_para<uni10_double64>*  ) = {init_nsy, init_bsy};
    static void (*init_z[])(U_para<uni10_complex128>*) = {init_nsy, init_bsy};

    static void (*init_val_d[])(uni10_double64,  U_para<uni10_double64>*  ) = {init_val_nsy, init_val_bsy};
    static void (*init_val_z[])(uni10_complex128,U_para<uni10_complex128>*) = {init_val_nsy, init_val_bsy};


    //static void (*grouping_d[])(U_para<uni10_double64>*  ) = {grouping_nsy, grouping_bsy};
    //static void (*grouping_z[])(U_para<uni10_complex128>*) = {grouping_nsy, grouping_bsy};

    //static void (*initBlocks_d[])(U_para<uni10_double64>*  ) = {initBlocks_nsy, initBlocks_bsy};
    //static void (*initBlocks_z[])(U_para<uni10_complex128>*) = {initBlocks_nsy, initBlocks_bsy};
    
    // link general variables to structure variables;
    template<typename U, typename T, typename... Args>
      void var_meta_links(U style, T* _para, Args... args){
        meta_link();
      }

    // Function overload for UniTensor<T>::init();
    void init(U_para<uni10_double64>*   para, contain_type _style);
    void init(U_para<uni10_complex128>* para, contain_type _style);

    // Function overload for UniTensor<T>::UniTensor(T val, contain_type style);
    void init_val(uni10_double64 val,   U_para<uni10_double64>*   para, contain_type _style);
    void init_val(uni10_complex128 val, U_para<uni10_complex128>* para, contain_type _style);


  };

};

#endif
