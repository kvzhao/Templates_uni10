#include "uni10/uni10_tools/uni10_tools_cpu.h"

namespace uni10{

  void* uni10_elem_alloc(uni10_int memsize){
    void* ptr = NULL;
    ptr = malloc(memsize);

    uni10_error_msg(ptr==NULL,"Fails in allocating memory.");

    return ptr;
  }

  void* uni10_elem_copy(void* des, const void* src, uni10_int memsize){
    return memcpy(des, src, memsize);
  }

  void uni10_elem_free(void* ptr, uni10_int memsize){
    free(ptr);
    ptr = NULL;
  }

  void uni10_elemBzero(void* ptr, uni10_int memsize){
    memset(ptr, 0, memsize);
  }

  // For double 
  //
  void uni10_setDiag(uni10_double* elem, uni10_double* diag_elem, uni10_int m, uni10_int n, uni10_int diag_n){

    uni10_int min = m < n ? m : n;

    min = min < diag_n ? min : diag_n;

    for(uni10_int i = 0; i < min; i++)
      elem[i * n + i] = diag_elem[i];

  }
  void uni10_getDiag(uni10_double* elem, uni10_double* diag_elem, uni10_int m, uni10_int n, uni10_int diag_n){

    uni10_int min = m < n ? m : n;

    min = min < diag_n ? min : diag_n;

    for(uni10_int i = 0; i < min; i++)
      diag_elem[i] = elem[i * n + i];

  }

  // For complex 
  //
  void uni10_setDiag(uni10_complex* elem, uni10_complex* diag_elem, uni10_int m, uni10_int n, uni10_int diag_n){

    uni10_int min = m < n ? m : n;

    min = min < diag_n ? min : diag_n;

    for(uni10_int i = 0; i < min; i++)
      elem[i * n + i] = diag_elem[i];

  }
  void uni10_getDiag(uni10_complex* elem, uni10_complex* diag_elem, uni10_int m, uni10_int n, uni10_int diag_n){

    uni10_int min = m < n ? m : n;

    min = min < diag_n ? min : diag_n;

    for(uni10_int i = 0; i < min; i++)
      diag_elem[i] = elem[i * n + i];

  }


} /* namespace uni10 */
