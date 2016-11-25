#include "uni10/uni10_api/tensor_tools.h"

namespace uni10{
  
  namespace tensor_tools{

    // Overload for UniTensor<T>::init_para();
    U_para<uni10_double64>* init_para(U_para<uni10_double64>* para, contain_type style){


      para = init_para_d[style](para);
      return para;

    }

    U_para<uni10_complex128>* init_para(U_para<uni10_complex128>* para, contain_type style){

      para = init_para_z[style](para);
      return para;

    }

    // Overload for UniTensor<T>::init_para();
    void free_para(U_para<uni10_double64>* para, contain_type style){

      para = init_para_d[style](para);

    }

    void free_para(U_para<uni10_complex128>* para, contain_type style){

      para = init_para_z[style](para);

    }

    // Overload for UniTensor<T>::init();
    void init(U_para<uni10_double64>* para, contain_type style){

      init_d[style](para);

    }

    void init(U_para<uni10_complex128>* para, contain_type style){

      init_z[style](para);

    }

  };
};
