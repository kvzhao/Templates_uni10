#ifndef __UNI10_TYPE_H__
#define __UNI10_TYPE_H__

#include <stdio.h>
#include <string.h>

#include <complex>

#define CPU 1
// The excute type of Uni10.
#ifdef CPU
#define _type _cpu
#elif  GPU
#define _type _gpu
#elif  MPI
#define _type _mpi
#endif


typedef int                    uni10_int32;       // shout int.
typedef long long int          uni10_int64;       // MKL uses long long int, not int64_t.
typedef size_t                 uni10_uint64;      
typedef float                  uni10_float32;
typedef double                 uni10_double64;
typedef std::complex<float>    uni10_complex64;
typedef std::complex<double>   uni10_complex128;

typedef int                    uni10_exu_type;    // To store the exu_type.

// Generate the typename of the uni10 system information.
#define info_type_helper(sysinfo, type)  sysinfo##type
#define info_type(sysinfo, type)  info_type_helper(sysinfo, type)

// Generate the suitable function name.
#define uni10_func_helper(func, type)  func##type
#define uni10_func(func, type)  uni10_func_helper(func, type)

// Generate the class name of uni10_elem.
#define elem_type_helper(UNI10ELEM, type)  UNI10ELEM##type
#define elem_type(UNI10ELEM, type)  elem_type_helper(UNI10ELEM, type)

//typedef double*               uni10_double_ptr;
//typedef std::complex<double>* uni10_complex_ptr;

#endif
