#include "uni10/uni10_env_info.h"
#include "uni10/uni10_lapack_cpu/tools_lapack_cpu/uni10_tools_cpu.h"

namespace uni10{

  void* uni10_elem_alloc_cpu(uni10_uint64 memsize){

    void* ptr = NULL;
    ptr = malloc(memsize);

    env_variables.use_memsize(memsize);
    uni10_error_msg(ptr==NULL, "%s","Fails in allocating memory.");

    return ptr;
  }

  void* uni10_elem_copy_cpu(void* des, const void* src, uni10_uint64 memsize){
    return memcpy(des, src, memsize);
  }

  void uni10_elem_free_cpu(void* ptr, uni10_uint64 memsize){

    free(ptr);
    env_variables.use_memsize(-memsize);
    ptr = NULL;

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

  void uni10_getUpTri_cpu(uni10_double64* elem, uni10_double64* tri_elem, uni10_uint64 m, uni10_uint64 n){

    uni10_uint64 min = m < n ? m : n;

    for(uni10_uint64 i = 0; i < min; i++)
      uni10_elem_copy_cpu(tri_elem + i*min+i, elem + i*n + (n-min)+i, (min-i)*sizeof(uni10_double64));

  }

  void uni10_getDnTri_cpu(uni10_double64* elem, uni10_double64* tri_elem, uni10_uint64 m, uni10_uint64 n){

    uni10_uint64 min = m < n ? m : n;

    for(uni10_uint64 i = 0; i < min; i++)
      uni10_elem_copy_cpu(tri_elem + i*min, elem + i*(m-min)+i*n, (i+1)*sizeof(uni10_double64));

  }

  void uni10_print_elem_i(const uni10_double64& elem_i){

    fprintf(stdout, " %8.4f", elem_i);

  }

  uni10_double64 UNI10_REAL( uni10_double64 elem_i ){
    return elem_i;
  }

  void uni10_getUpTri_cpu(uni10_complex128* elem, uni10_complex128* tri_elem, uni10_uint64 m, uni10_uint64 n){

    uni10_error_msg(m < n, "%s", "Can't get the upper triangular matrix when Rnum < Cnum.");
    uni10_uint64 min = m < n ? m : n;

    for(uni10_uint64 i = 0; i < min; i++)
      uni10_elem_copy_cpu(tri_elem + i*min+i, elem + i*n+(n-min)+i, (min-i)*sizeof(uni10_complex128));

  }

  void uni10_getDnTri_cpu(uni10_complex128* elem, uni10_complex128* tri_elem, uni10_uint64 m, uni10_uint64 n){

    uni10_error_msg(m > n, "%s", "Can't get the lower triangular matrix when Rnum < Cnum.");
    uni10_uint64 min = m < n ? m : n;

    for(uni10_uint64 i = 0; i < min; i++)
      uni10_elem_copy_cpu(tri_elem + i*min, elem + i*(m-min)+i*n, (i+1)*sizeof(uni10_complex128));

  }

  uni10_double64 UNI10_IMAG( uni10_double64 elem_i ){
    return 0.;
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

  void uni10_print_elem_i(const uni10_complex128& elem_i){

    fprintf(stdout, " %8.4f+%8.4fi", Z_REAL( elem_i ), Z_IMAG( elem_i ) );

  }

  uni10_double64 UNI10_REAL( uni10_complex128 elem_i ){
    return elem_i.real();
  }

  uni10_double64 UNI10_IMAG( uni10_complex128 elem_i ){
    return elem_i.imag();
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

  void shrinkWithoutFree(uni10_uint64 memsize){

    env_variables.use_memsize(-memsize);

  }

} /* namespace uni10 */
