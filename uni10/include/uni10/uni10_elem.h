#ifndef __UNI10_ELEM_H__
#define __UNI10_ELEM_H__

//#define CPU
//#define LAPACK

#if defined(LAPACK) && defined(CPU)
#include "uni10/uni10_lapack_cpu/uni10_elem_lapack_cpu.h"
#endif

#if defined(MAGMA) && defined(GPU)
#include "uni10/uni10_magma_cpu/uni10_elem_magma_gpu.h"
#endif

#endif
