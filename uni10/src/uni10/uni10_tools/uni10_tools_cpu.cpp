#include "uni10/uni10_env_info.h"
#include "uni10/uni10_tools/uni10_tools_cpu.h"

namespace uni10{

  void* uni10_elem_alloc_cpu(uni10_uint64 memsize){

    void* ptr = NULL;
    ptr = malloc(memsize);

    env_variables.use_memsize(memsize);
    uni10_error_msg(ptr==NULL,"Fails in allocating memory.");

    return ptr;
  }

  void* uni10_elem_copy_cpu(void* des, const void* src, uni10_uint64 memsize){
    return memcpy(des, src, memsize);
  }

  void uni10_elem_free_cpu(void* ptr, uni10_uint64 memsize){

    free(ptr);
    ptr = NULL;
    env_variables.use_memsize(-memsize);

  }

  void uni10_elemBzero_cpu(void* ptr, uni10_uint64 memsize){
    memset(ptr, 0, memsize);
  }

  // For double 
  //
  void uni10_setDiag_cpu(uni10_double64* elem, uni10_double64* diag_elem, uni10_uint64 m, uni10_uint64 n, uni10_uint64 diag_n){

    uni10_uint64 min = m < n ? m : n;

    min = min < diag_n ? min : diag_n;

    for(uni10_uint64 i = 0; i < min; i++)
      elem[i * n + i] = diag_elem[i];

  }
  void uni10_getDiag_cpu(uni10_double64* elem, uni10_double64* diag_elem, uni10_uint64 m, uni10_uint64 n, uni10_uint64 diag_n){

    uni10_uint64 min = m < n ? m : n;

    min = min < diag_n ? min : diag_n;

    for(uni10_uint64 i = 0; i < min; i++)
      diag_elem[i] = elem[i * n + i];

  }

  // For complex 
  //
  void uni10_setDiag_cpu(uni10_complex128* elem, uni10_complex128* diag_elem, uni10_uint64 m, uni10_uint64 n, uni10_uint64 diag_n){

    uni10_uint64 min = m < n ? m : n;

    min = min < diag_n ? min : diag_n;

    for(uni10_uint64 i = 0; i < min; i++)
      elem[i * n + i] = diag_elem[i];

  }
  void uni10_getDiag_cpu(uni10_complex128* elem, uni10_complex128* diag_elem, uni10_uint64 m, uni10_uint64 n, uni10_uint64 diag_n){

    uni10_uint64 min = m < n ? m : n;

    min = min < diag_n ? min : diag_n;

    for(uni10_uint64 i = 0; i < min; i++)
      diag_elem[i] = elem[i * n + i];

  }
  
  // Convert
  void uni10_elem_cast_cpu(uni10_complex128* des, uni10_double64* src, uni10_uint64 N){

    for(uni10_uint64 i = 0; i < N; i++)
      des[i] = src[i];

  }

  void uni10_elem_cast_cpu(uni10_double64* des, uni10_complex128* src, uni10_uint64 N){

    for(uni10_uint64 i = 0; i < N; i++)
      des[i] = src[i].real();

  }

} /* namespace uni10 */
