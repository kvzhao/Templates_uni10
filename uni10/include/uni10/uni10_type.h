#ifndef __UNI10_TYPE_H__
#define __UNI10_TYPE_H__

#include <stdio.h>
#include <string.h>

#include <complex>

typedef long long int          uni10_int64;       // MKL uses long long int, not int64_t.
typedef size_t                 uni10_uint64;      
typedef float                  uni10_float32;
typedef double                 uni10_double64;
typedef std::complex<float>    uni10_complex64;
typedef std::complex<double>   uni10_complex128;

typedef int                    uni10_exu_type;    // To store the exu_type.

//typedef double*               uni10_double_ptr;
//typedef std::complex<double>* uni10_complex_ptr;

#endif
