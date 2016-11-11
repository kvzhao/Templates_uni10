#ifndef __UNI10_ELEM_LINALG_H__
#define __UNI10_ELEM_LINALG_H__

//#define CPU
//#define LAPACK

#if defined(CPU) && defined(LAPACK)
#include "uni10/uni10_lapack_cpu/uni10_elem_linalg_lapack_cpu.h"
#endif

#endif

