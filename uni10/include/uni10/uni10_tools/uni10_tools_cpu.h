#ifndef __UNI10_TOOLS_CPU_H__
#define __UNI10_TOOLS_CPU_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "uni10/uni10_type.h"
#include "uni10/uni10_error.h"

namespace uni10{

  void* uni10_elem_alloc(uni10_int memsize);

  void* uni10_elem_copy(void* des, const void* src, uni10_int memsize);

  void uni10_elem_free(void* ptr, uni10_int memsize);

  void uni10_elemBzero(void* ptr, uni10_int memsize);

  //For double ptr.
  //
  void uni10_setDiag(uni10_double* elem, uni10_double* diag_elem, uni10_int M, uni10_int N, uni10_int diag_N);

  void uni10_getDiag(uni10_double* elem, uni10_double* diag_elem, uni10_int M, uni10_int N, uni10_int diag_N);

  //For complex ptr.
  //
  void uni10_setDiag(uni10_complex* elem, uni10_complex* diag_elem, uni10_int M, uni10_int N, uni10_int diag_N);

  void uni10_getDiag(uni10_complex* elem, uni10_complex* diag_elem, uni10_int M, uni10_int N, uni10_int diag_N);

}/* namespace uni10 */
#endif /* UNI10_TOOLS_H */
