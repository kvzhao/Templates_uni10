#include "uni10/uni10_api/tensor_tools.h"

namespace uni10{
  
  namespace tensor_tools{

    // Overload for UniTensor<T>::init();
    void init(U_para<uni10_double64>* para, contain_type style){

      init_d[style](para);

    }

    void init(U_para<uni10_complex128>* para, contain_type style){

      init_z[style](para);

    }

    // Overload for UniTensor<T>::UniTensor(T val, contain_type style);
    void init_val(uni10_double64 val, U_para<uni10_double64>* para, contain_type style){

      init_val_d[style](val, para);

    }

    void init_val(uni10_complex128 val, U_para<uni10_complex128>* para, contain_type style){

      init_val_z[style](val, para);

    }

  };
};
