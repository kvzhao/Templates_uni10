#ifndef __UNI10_TENSOR_TOOLS_H__
#define __UNI10_TENSOR_TOOLS_H__

#include <stdio.h>

#include "uni10/uni10_api/nsy_tensor_tools.h"
#include "uni10/uni10_api/bsy_tensor_tools.h"

namespace uni10{

  namespace tensor_tools{
    // The function tools to initialize UniTensor.
    static U_para<uni10_double64>*  (*init_para_d[])(U_para<uni10_double64>*  ) = {init_para_nsy, init_para_bsy};
    static U_para<uni10_complex128>*(*init_para_z[])(U_para<uni10_complex128>*) = {init_para_nsy, init_para_bsy};

    // The function tools to initialize UniTensor.
    static void (*free_para_d[])(U_para<uni10_double64>*  ) = {free_para_nsy, free_para_bsy};
    static void (*free_para_z[])(U_para<uni10_complex128>*) = {free_para_nsy, free_para_bsy};

    static void (*init_d[])(U_para<uni10_double64>*  ) = {init_nsy, init_bsy};
    static void (*init_z[])(U_para<uni10_complex128>*) = {init_nsy, init_bsy};

    //static void (*grouping_d[])(U_para<uni10_double64>*  ) = {grouping_nsy, grouping_bsy};
    //static void (*grouping_z[])(U_para<uni10_complex128>*) = {grouping_nsy, grouping_bsy};

    //static void (*initBlocks_d[])(U_para<uni10_double64>*  ) = {initBlocks_nsy, initBlocks_bsy};
    //static void (*initBlocks_z[])(U_para<uni10_complex128>*) = {initBlocks_nsy, initBlocks_bsy};
    
    // Function overload for initialize U_paras;
    U_para<uni10_double64>*   init_para(U_para<uni10_double64>*   para, contain_type _style);
    U_para<uni10_complex128>* init_para(U_para<uni10_complex128>* para, contain_type _style);

    // Function overload for free parameters;
    void free_para(U_para<uni10_double64>*   para, contain_type _style);
    void free_para(U_para<uni10_complex128>* para, contain_type _style);

    // Function overload for UniTensor<T>::init();
    void init(U_para<uni10_double64>*   para, contain_type _style);
    void init(U_para<uni10_complex128>* para, contain_type _style);

  };

};

#endif
